use crate::hydro::Eos;

use super::{time::schemes::Scheme, Constraint};

pub const DIM: usize = 3;
pub type Arr<const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    [[[[f64; F]; VX]; VY]; VZ];
pub type BArr<const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    Box<Arr<F, VX, VY, VZ>>;

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint,
}

pub type Boundary<'a, const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    &'a (dyn Fn([i32; DIM], &Arr<F, VX, VY, VZ>) -> [f64; F] + Sync);
pub type Fun<
    'a,
    Opt,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> = &'a (dyn Fn(
    &Arr<F, VX, VY, VZ>,
    [&Arr<F, VX, VY, VZ>; 2],
    [&Arr<C, VX, VY, VZ>; 2],
    Constraint<F, C>,
    Boundary<F, VX, VY, VZ>,
    [i32; DIM], // position in index [x,y]
    [f64; DIM], // [dx, dy, dz]
    [f64; 2],   // [old t, current t]
    [f64; 2],   // [dt, current dt]
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
    const VZ: usize,
    const S: usize,
> {
    pub fun: Fun<'a, Opt, F, C, VX, VY, VZ>,
    pub constraints: Constraint<'a, F, C>,
    pub boundary: Boundary<'a, F, VX, VY, VZ>,
    pub post_constraints: Option<Constraint<'a, F, C>>,
    pub local_interaction: [i32; DIM],
    pub vstrs: (BArr<F, VX, VY, VZ>, BArr<C, VX, VY, VZ>),
    pub ovstrs: (BArr<F, VX, VY, VZ>, BArr<C, VX, VY, VZ>), // old
    pub total_diff_vs: BArr<F, VX, VY, VZ>,
    pub k: Box<[Arr<F, VX, VY, VZ>; S]>,
    pub r: Scheme<S>,
    pub dt: f64,
    pub dxs: [f64; DIM],
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
