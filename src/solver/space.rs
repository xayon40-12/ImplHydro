pub mod kl;
pub mod kt;

use crate::solver::{context::Boundary, Transform};

use self::{kl::kl, kt::kt};

pub type Flux<'a, const C: usize> = &'a dyn Fn(f64, [f64; C]) -> f64;
#[derive(Clone, Copy)]
pub enum Eigenvalues<'a, const C: usize> {
    Analytical(&'a (dyn Fn(f64, [f64; C]) -> f64 + Sync)), // Fn(time, non dymical variables) -> max eigenvalue
    ApproxConstraint(&'a (dyn Fn(f64, [f64; C], f64) -> f64 + Sync)), // Fn(time, non dynamical variables, approx eigenvalue) -> constrained eigenvalue
}

pub enum Dir {
    X,
    Y,
}

pub fn id_flux_limiter<const F: usize>(_t: f64, v: [f64; F]) -> [f64; F] {
    v
}

pub type SpaceDiff<'a, const F: usize, const VX: usize, const VY: usize, const C: usize> =
    &'a dyn Fn(
        &[[[f64; F]; VX]; VY],
        &[Boundary; 2],
        [i32; 2],
        Dir,
        f64,
        [Flux<C>; F],
        Transform<F, F>,
        Transform<F, C>,
        Eigenvalues<C>,
        Transform<F, F>,
        Transform<F, F>,
        f64,
        f64,
    ) -> [f64; F];

#[derive(Debug, Clone, Copy)]
pub enum Order {
    Order2,
    Order2Cut(f64),
    Order3,
    Order3Cut(f64),
}

pub fn order<'a, const F: usize, const VX: usize, const VY: usize, const C: usize>(
    o: Order,
) -> SpaceDiff<'a, F, VX, VY, C> {
    match o {
        Order::Order2Cut(_) | Order::Order2 => &kt,
        Order::Order3Cut(_) | Order::Order3 => &kl,
    }
}
