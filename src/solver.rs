pub mod context;
pub mod explicit;
pub mod fixpoint;
pub mod kt;
pub mod newton;
pub mod schemes;
pub mod utils;

use {
    context::{Context, Integration},
    explicit::explicit,
    fixpoint::fixpoint,
};

pub type Constraints<'a, const F: usize, const C: usize> =
    &'a (dyn Fn([f64; F]) -> [f64; C] + Sync);

pub fn save<const F: usize, const C: usize, const VX: usize, const VY: usize>(
    v: &[[[f64; F]; VX]; VY],
    constraints: Constraints<F, C>,
    names: &[&str; C],
    name: &str, // simulation name
    elapsed: f64,
    t0: f64,
    tend: f64,
    t: f64,
    dx: f64,
    maxdt: f64,
    cost: usize,
    integration: Integration,
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
            let mut s = format!("{:e} {:e}", x, y);
            for c in 0..C {
                s = format!("{} {:e}", s, vars[c]);
            }
            s = format!("{}\n", s);
            res = format!("{}{}", res, s);
        }
        res = format!("{}\n", res);
    }

    let dir = &format!("results/{}/{:e}", name, t);
    std::fs::create_dir_all(dir)?;
    std::fs::write(&format!("{}/data.txt", dir), res.as_bytes())?;
    let info = format!(
        "elapsed: {:e}\nt0: {:e}\ntend: {:e}\nt: {:e}\ncost: {}\nnx: {}\nny: {}\ndx: {:e}\nmaxdt: {:e}\nintegration: {:?}",
        elapsed, t0, tend, t, cost, VX, VY, dx, maxdt, integration
    );
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
    mut context: Context<Opt, F, C, VX, VY, S>,
    name: &str,
    integration: Integration,
    names: &[&str; C],
    constraints: Constraints<F, C>,
) -> ([[[f64; F]; VX]; VY], f64, usize, usize) {
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let now = Instant::now();

    let next_save = context.tend;
    let mut ctx_dt = None;
    while context.t < context.tend {
        let d = next_save - context.t;
        if d < 0.0 || d > 2.0 * context.dt {
            if let Some(dt) = ctx_dt {
                context.dt = dt;
                ctx_dt = None;
            }
        } else if d > context.dt {
            ctx_dt = Some(context.dt);
            context.dt = d / 2.0;
        } else {
            ctx_dt = Some(context.dt);
            context.dt = d;
        }
        tsteps += 1;
        let c = match integration {
            Integration::Explicit => explicit(&mut context),
            Integration::FixPoint => fixpoint(&mut context),
        };

        cost += c;
    }
    let cost = cost as usize;
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Elapsed: {:.2?}", elapsed);
    let err = save(
        &context.vs,
        constraints,
        names,
        name,
        elapsed,
        context.t0,
        context.tend,
        context.t,
        context.dx,
        context.maxdt,
        cost,
        integration,
    );
    let _elapsed = now.elapsed().as_secs();
    eprintln!("Elapsed with save: {:.2?}", _elapsed);
    match err {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    (context.vs, context.t, cost, tsteps)
}
