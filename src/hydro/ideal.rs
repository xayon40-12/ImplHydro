use crate::hydro::{Eos, VOID};

pub mod ideal1d;
pub mod ideal2d;

pub fn init_from_energy_1d<'a, const VX: usize>(
    t0: f64,
    es: [f64; VX],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn(usize, f64) -> [f64; 2] + 'a> {
    Box::new(move |i, _| {
        let e = es[i].max(VOID);
        let vars = [e, p(e), dpde(e), 1.0, 0.0];
        ideal1d::f0(t0, vars)
    })
}

pub fn init_from_energy_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    es: [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; 3] + 'a> {
    Box::new(move |(i, j), _| {
        let e = es[j][i].max(VOID);
        let vars = [e, p(e), dpde(e), 1.0, 0.0, 0.0];
        ideal2d::f0(t0, vars)
    })
}
