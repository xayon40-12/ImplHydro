use crate::hydro::Eos;

use super::{time::schemes::Scheme, Constraint};

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint(Option<f64>), // takes a possible constant to multiply to dt^p for error checking in the iterative solver
}

pub type Boundary<'a> = &'a (dyn Fn(i32, usize) -> usize + Sync);
pub type Fun<'a, Opt, const F: usize, const C: usize, const VX: usize, const VY: usize> =
    &'a (dyn Fn(
        [&[[[f64; F]; VX]; VY]; 2],
        [&[[[f64; C]; VX]; VY]; 2],
        Constraint<F, C>,
        &[Boundary; 2],
        [i32; 2], // position in index [x,y]
        f64,      // dx
        [f64; 2], // [old t, current t]
        [f64; 2], // [dt, current dt]
        &Opt,
    ) -> [f64; F]
             + Sync);

#[derive(Clone)]
pub struct Context<
    'a,
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
> {
    pub fun: Fun<'a, Opt, F, C, VX, VY>,
    pub constraints: Constraint<'a, F, C>,
    pub boundary: &'a [Boundary<'a>; 2],
    pub post_constraints: Option<Constraint<'a, F, C>>,
    pub local_interaction: [i32; 2],
    pub vstrs: ([[[f64; F]; VX]; VY], [[[f64; C]; VX]; VY]),
    pub ovstrs: ([[[f64; F]; VX]; VY], [[[f64; C]; VX]; VY]), // old
    pub total_diff_vs: [[[f64; F]; VX]; VY],
    pub k: [[[[f64; F]; VX]; VY]; S],
    pub r: Scheme<S>,
    pub dt: f64,
    pub dx: f64,
    pub maxdt: f64,
    pub er: f64,
    pub t: f64,
    pub ot: f64,
    pub t0: f64,
    pub tend: f64,
    pub opt: Opt,
    pub p: Eos<'a>,
    pub dpde: Eos<'a>,
    pub freezeout_energy: Option<f64>,
}
