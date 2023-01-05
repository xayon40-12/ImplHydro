use crate::hydro::Pressure;

use super::{schemes::Scheme, Constraints};

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    Explicit,
    FixPoint,
}

#[derive(Debug, Clone, Copy)]
pub enum ToCompute {
    Integrated,
    NonIntegrated,
    All,
}

impl ToCompute {
    pub fn integrated(&self) -> bool {
        match self {
            ToCompute::Integrated => true,
            _ => false,
        }
    }
    pub fn nonintegrated(&self) -> bool {
        match self {
            ToCompute::NonIntegrated => true,
            _ => false,
        }
    }
    pub fn all(&self) -> bool {
        match self {
            ToCompute::All => true,
            _ => false,
        }
    }
}

pub type Boundary<'a> = &'a (dyn Fn(i32, usize) -> usize + Sync);
pub type Fun<'a, Opt, const F: usize, const C: usize, const VX: usize, const VY: usize> =
    &'a (dyn Fn(
        [&[[[f64; F]; VX]; VY]; 2],
        Constraints<F, C>,
        &[Boundary; 2],
        [i32; 2], // position in index
        f64,      // dx
        f64,      // er
        [f64; 2], // [old t, current t]
        [f64; 2], // [dt, current dt]
        &Opt,
        ToCompute,
        Pressure,
        Pressure,
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
    pub constraints: Constraints<'a, F, C>,
    pub boundary: &'a [Boundary<'a>; 2],
    pub local_interaction: [i32; 2],
    pub vs: [[[f64; F]; VX]; VY],
    pub k: [[[[f64; F]; VX]; VY]; S],
    pub integrated: [bool; F],
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
