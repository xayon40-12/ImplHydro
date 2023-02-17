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
        ];
        [
            shear2d::fi00(t0, vars),
            shear2d::fi01(t0, vars),
            shear2d::fi02(t0, vars),
            shear2d::u0pi00(t0, vars),
            shear2d::u0pi01(t0, vars),
            shear2d::u0pi02(t0, vars),
            shear2d::u0pi11(t0, vars),
            shear2d::u0pi12(t0, vars),
            shear2d::u0pi22(t0, vars),
        ]
    })
}
