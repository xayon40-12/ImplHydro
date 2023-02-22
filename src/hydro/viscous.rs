use crate::hydro::{Eos, VOID};

use super::F_SHEAR_2D;

pub mod shear2d;

pub fn init_from_energy_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    es: [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_SHEAR_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let e = es[j][i].max(VOID);
        let vars = [
            e,
            p(e),
            dpde(e),
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        shear2d::fitutpi(t0, vars)
    })
}
