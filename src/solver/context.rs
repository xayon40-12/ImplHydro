use crate::{hydro::Eos, FLOAT};

use super::{time::schemes::Scheme, Constraint};

pub const DIM: usize = 3;
pub type Arr<const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    [[[[FLOAT; F]; VX]; VY]; VZ];
pub type BArr<const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    Box<Arr<F, VX, VY, VZ>>;

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint,
}

pub type Boundary<'a, const F: usize, const VX: usize, const VY: usize, const VZ: usize> =
    &'a (dyn Fn([i32; DIM], &Arr<F, VX, VY, VZ>) -> [FLOAT; F] + Sync);
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
    [i32; DIM],   // position in index [x,y]
    [FLOAT; DIM], // [dx, dy, dz]
    [FLOAT; 2],   // [old t, current t]
    [FLOAT; 2],   // [dt, current dt]
    &Opt,
) -> [FLOAT; F]
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
    pub dt: FLOAT,
    pub dxs: [FLOAT; DIM],
    pub maxdt: FLOAT,
    pub t: FLOAT,
    pub ot: FLOAT,
    pub t0: FLOAT,
    pub tend: FLOAT,
    pub opt: Opt,
    pub p: Eos<'a>,
    pub dpde: Eos<'a>,
    pub freezeout_energy: Option<FLOAT>,
}
