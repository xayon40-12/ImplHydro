use crate::newton::{newton, Context};

pub fn run<const F: usize, const VX: usize, const VY: usize, const S: usize>(
    mut context: Context<F, VX, VY, S>,
) -> ([[[f64; VX]; VY]; F], usize, usize) {
    let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let tsteps = ((context.tend - context.tbeg) / context.dt) as usize;
    let mut tot_f = 0;
    for _i in 0..tsteps {
        let c = newton(&mut context);
        let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        tot_f += c;
        // if _i % (1 + tsteps / 10) == 0 {
        //     println!("i: {}, c: {}", _i, c);
        // }
    }
    (context.vs, tot_f, tsteps)
}
