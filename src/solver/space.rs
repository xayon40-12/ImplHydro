pub mod kt;

pub type Flux<'a, const F: usize, const C: usize> = &'a dyn Fn(f64, [f64; C]) -> [f64; F];
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
