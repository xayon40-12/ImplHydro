use crate::{fixedpoint::fixedpoint, newton::Context};

pub fn run<Opt: Sync, const F: usize, const VX: usize, const VY: usize, const S: usize>(
    mut context: Context<Opt, F, VX, VY, S>,
) -> ([[[f64; F]; VX]; VY], usize, usize) {
    let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut tot_f = 0;
    let mut tsteps = 0;
    let tend = context.tend;
    let mut t = context.t;
    while t < tend {
        tsteps += 1;
        let c = fixedpoint(&mut context);
        t = context.t;
        let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        tot_f += c;
        // if _i % (1 + tsteps / 10) == 0 {
        //     println!("i: {}, c: {}", _i, c);
        // }
    }
    (context.vs, tot_f, tsteps)
}
