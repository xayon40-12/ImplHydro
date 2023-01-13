pub mod kl;
pub mod kt;

use crate::solver::{context::Boundary, Transform};

use self::{kl::kl, kt::kt};

pub type Flux<'a, const C: usize> = &'a dyn Fn(f64, [f64; C]) -> f64;
pub type Eigenvalues<'a, const C: usize> = &'a dyn Fn(f64, [f64; C]) -> f64;

pub enum Dir {
    X,
    Y,
}

pub fn id_flux_limiter<const F: usize>(_t: f64, v: [f64; F]) -> [f64; F] {
    v
}

pub type SpaceDiff<
    'a,
    const F: usize,
    const VX: usize,
    const VY: usize,
    const C: usize,
    const N: usize,
> = &'a dyn Fn(
    &[[[f64; F]; VX]; VY],
    &[Boundary; 2],
    [i32; 2],
    Dir,
    f64,
    [Flux<C>; N],
    Transform<F, F>,
    Transform<F, C>,
    Eigenvalues<C>,
    Transform<F, F>,
    Transform<F, F>,
    f64,
    f64,
) -> [f64; N];

pub enum Order {
    O2,
    O3,
}

pub fn order<
    'a,
    const F: usize,
    const VX: usize,
    const VY: usize,
    const C: usize,
    const N: usize,
>(
    o: Order,
) -> SpaceDiff<'a, F, VX, VY, C, N> {
    match o {
        Order::O2 => &kt,
        Order::O3 => &kl,
    }
}
