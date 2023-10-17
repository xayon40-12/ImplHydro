use super::{
    context::{Arr, Boundary, DIM},
    Constraint, Transform,
};

pub mod kt;

pub type SpaceDerivative<
    'a,
    const N: usize,
    const F: usize,
    const C: usize,
    const SD: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> = &'a dyn Fn(
    (&Arr<F, VX, VY, VZ>, &Arr<C, VX, VY, VZ>), // (vrs,trs)
    Boundary<F, VX, VY, VZ>,                    // bound
    [i32; DIM],                                 // [x,y,z]
    f64,                                        // t
    FluxInfos<N, F, C, SD>,                     // flux_infos
    Constraint<F, C>,                           // constraits
    Transform<F>,                               // pre_flux_limiter
    Transform<F>,                               // post_flux_limiter
    f64,                                        // dx
    f64,                                        // theta
) -> [InDir<([f64; F], [f64; SD])>; N];
pub type Flux<'a, const F: usize, const C: usize> = &'a dyn Fn(f64, [f64; C]) -> [f64; F];
#[derive(Clone, Copy)]
pub enum Eigenvalues<'a, const C: usize> {
    Analytical(&'a (dyn Fn(f64, [f64; C]) -> f64 + Sync)), // Fn(time, non dymical variables) -> max eigenvalue
    ApproxConstraint(&'a (dyn Fn(f64, [f64; C], f64) -> f64 + Sync)), // Fn(time, non dynamical variables, approx eigenvalue) -> constrained eigenvalue
}

pub type FluxInfos<'a, const N: usize, const F: usize, const C: usize, const SD: usize> =
    [InDir<FluxInfo<'a, F, C, SD>>; N];

pub struct FluxInfo<'a, const F: usize, const C: usize, const SD: usize> {
    pub flux: Flux<'a, F, C>,
    pub secondary: Flux<'a, SD, C>,
    pub eigenvalues: Eigenvalues<'a, C>,
}

pub enum Dir {
    X,
    Y,
    Z,
}

#[derive(Copy, Clone)]
pub enum InDir<T> {
    X(T),
    Y(T),
    Z(T),
}

impl<T> InDir<T> {
    pub fn dir(&self) -> Dir {
        match self {
            InDir::X(_) => Dir::X,
            InDir::Y(_) => Dir::Y,
            InDir::Z(_) => Dir::Z,
        }
    }
    pub fn iner(self) -> T {
        match self {
            InDir::X(a) => a,
            InDir::Y(a) => a,
            InDir::Z(a) => a,
        }
    }
    pub fn map<U, F: Fn(T) -> U>(self, f: F) -> InDir<U> {
        match self {
            InDir::X(a) => InDir::X(f(a)),
            InDir::Y(a) => InDir::Y(f(a)),
            InDir::Z(a) => InDir::Z(f(a)),
        }
    }
}

pub fn id_flux_limiter<const F: usize>(_t: f64, v: [f64; F]) -> [f64; F] {
    v
}
