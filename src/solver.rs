use crate::{context::Context, explicit::explicit, fixedpoint::fixedpoint};

pub type Constraints<'a, const F: usize, const C: usize> = &'a dyn Fn([f64; F]) -> [f64; C];

pub fn save<const F: usize, const C: usize, const VX: usize, const VY: usize>(
    v: &[[[f64; F]; VX]; VY],
    constraints: Constraints<F, C>,
    names: &[&str; C],
    t: f64,
    dx: f64,
    cost: usize,
) -> std::io::Result<()> {
    let mut res = format!("# t {:e}\n# cost {}\n# x y", t, cost);
    for c in 0..C {
        res = format!("{} {}", res, names[c]);
    }
    res = format!("{}\n", res);

    for j in 0..VY {
        for i in 0..VX {
            let vars = constraints(v[j][i]);
            let y = (j as f64 - ((VY - 1) as f64) / 2.0) * dx;
            let x = (i as f64 - ((VX - 1) as f64) / 2.0) * dx;
            let mut s = format!("{} {}", x, y);
            for c in 0..C {
                s = format!("{} {}", s, vars[c]);
            }
            s = format!("{}\n", s);
            res = format!("{}{}", res, s);
        }
        res = format!("{}\n", res);
    }

    let time = format!("{:e}", t);
    let dir = &format!("results/{}", time);
    std::fs::create_dir_all(dir)?;
    std::fs::write(&format!("{}/data.txt", dir), res.as_bytes())?;
    let info = format!("time: {}\ncost: {}\n", time, cost);
    std::fs::write(&format!("{}/info.txt", dir), info.as_bytes())?;

    Ok(())
}

pub fn run<
    Opt: Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    mut context: Context<Opt, F, VX, VY, S>,
    names: &[&str; C],
    constraints: Constraints<F, C>,
) -> ([[[f64; F]; VX]; VY], f64, usize, usize) {
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let now = Instant::now();
    while context.t < context.tend {
        tsteps += 1;
        let c = explicit(&mut context);
        cost += c;
    }
    let cost = cost as usize;
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    let err = save(&context.vs, constraints, names, context.t, context.dx, cost);
    match err {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    (context.vs, context.t, cost, tsteps)
}
