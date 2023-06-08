use crate::hydro::Eos;

use super::{time::schemes::Scheme, Constraint};

pub type Arr<const F: usize, const VX: usize, const VY: usize> = [[[f64; F]; VX]; VY];
pub type BArr<const F: usize, const VX: usize, const VY: usize> = Box<Arr<F, VX, VY>>;

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint,
}

pub type Boundary<'a, const F: usize, const VX: usize, const VY: usize> =
    &'a (dyn Fn([i32; 2], &Arr<F, VX, VY>) -> [f64; F] + Sync);
pub type Fun<'a, Opt, const F: usize, const C: usize, const VX: usize, const VY: usize> =
    &'a (dyn Fn(
        [&Arr<F, VX, VY>; 2],
        [&Arr<C, VX, VY>; 2],
        Constraint<F, C>,
        Boundary<F, VX, VY>,
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
    pub boundary: Boundary<'a, F, VX, VY>,
    pub post_constraints: Option<Constraint<'a, F, C>>,
    pub local_interaction: [i32; 2],
    pub vstrs: (Box<Arr<F, VX, VY>>, Box<Arr<C, VX, VY>>),
    pub ovstrs: (Box<Arr<F, VX, VY>>, Box<Arr<C, VX, VY>>), // old
    pub total_diff_vs: Box<Arr<F, VX, VY>>,
    pub k: Box<[Arr<F, VX, VY>; S]>,
    pub r: Scheme<S>,
    pub dt: f64,
    pub dx: f64,
    pub maxdt: f64,
    pub t: f64,
    pub ot: f64,
    pub t0: f64,
    pub tend: f64,
    pub opt: Opt,
    pub p: Eos<'a>,
    pub dpde: Eos<'a>,
    pub freezeout_energy: Option<f64>,
}
