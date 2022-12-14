use crate::newton::{newton, Context};

pub fn run<const F: usize, const V: usize, const S: usize>(
    mut context: Context<F, V, S>,
) -> ([[f64; V]; F], usize, usize) {
    let mul = context.local_interaction * 2 + 1;
    let tsteps = ((context.tend - context.tbeg) / context.dt) as usize;
    let mut tot_f = 0;
    for _i in 0..tsteps {
        let (nvs, c) = newton(&context);
        let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        context.vs = nvs;
        tot_f += c;
        // if _i % (1 + tsteps / 10) == 0 {
        //     println!("i: {}, c: {}", _i, c);
        // }
    }
    (context.vs, tot_f, tsteps)
}
