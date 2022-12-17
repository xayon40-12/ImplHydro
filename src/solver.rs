use crate::{fixedpoint::fixedpoint, newton::Context};

pub type Constraints<'a, const F: usize, const C: usize> = &'a dyn Fn([f64; F]) -> [f64; C];

pub fn save<const F: usize, const C: usize, const VX: usize, const VY: usize>(
    v: &[[[f64; F]; VX]; VY],
    constraints: Constraints<F, C>,
    names: &[&str; C],
    t: f64,
    dx: f64,
    iterations: usize,
) {
    let mut start = format!("# t {}\n# iterations {}\n# x y", t, iterations);
    for c in 0..C {
        start = format!("{} {}", start, names[c]);
    }
    println!("{}", start);

    for j in 0..VY {
        for i in 0..VX {
            let vars = constraints(v[j][i]);
            let y = (j as f64 - ((VY - 1) as f64) / 2.0) * dx;
            let x = (i as f64 - ((VX - 1) as f64) / 2.0) * dx;
            let mut s = format!("{} {}", x, y);
            for c in 0..C {
                s = format!("{} {}", s, vars[c]);
            }
            println!("{}", s);
        }
        println!("");
    }
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
) -> ([[[f64; F]; VX]; VY], usize, usize) {
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut tot_f = 0;
    let mut tsteps = 0;
    let tend = context.tend;
    let mut t = context.t;
    while t < tend {
        tsteps += 1;
        let c = fixedpoint(&mut context);
        t = context.t;
        // let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        tot_f += c;
        // if _i % (1 + tsteps / 10) == 0 {
        //     println!("i: {}, c: {}", _i, c);
        // }
    }
    save(
        &context.vs,
        constraints,
        names,
        context.t,
        context.dx,
        tot_f,
    );
    (context.vs, tot_f, tsteps)
}
