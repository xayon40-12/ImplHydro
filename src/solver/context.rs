use crate::hydro::Pressure;

use super::{schemes::Scheme, Transform};

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint,
}

pub type Boundary<'a> = &'a (dyn Fn(i32, usize) -> usize + Sync);
pub type Fun<'a, Opt, const F: usize, const C: usize, const VX: usize, const VY: usize> =
    &'a (dyn Fn(
        [&[[[f64; F]; VX]; VY]; 2],
        Transform<F, C>,
        &[Boundary; 2],
        [i32; 2], // position in index
        f64,      // dx
        [f64; 2], // [old t, current t]
        [f64; 2], // [dt, current dt]
        &Opt,
    ) -> [f64; F]
             + Sync);

pub struct Context<
    'a,
    Opt: Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
> {
    pub fun: Fun<'a, Opt, F, C, VX, VY>,
    pub constraints: Transform<'a, F, F>,
    pub transform: Transform<'a, F, C>,
    pub boundary: &'a [Boundary<'a>; 2],
    pub local_interaction: [i32; 2],
    pub vs: [[[f64; F]; VX]; VY],
    pub k: [[[[f64; F]; VX]; VY]; S],
    pub r: Scheme<S>,
    pub dt: f64,
    pub dx: f64,
    pub maxdt: f64,
    pub er: f64,
    pub t: f64,
    pub t0: f64,
    pub tend: f64,
    pub opt: Opt,
    pub p: Pressure<'a>,
    pub dpde: Pressure<'a>,
}
