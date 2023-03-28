use crate::{
    hydro::{Eos, VOID},
    solver::time::newton::newton,
};

use super::F_SHEAR_2D;

pub mod shear2d;

pub fn init_from_entropy_density_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    s: [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
    temperature: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_SHEAR_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let normalization = 50.0; // TODO adapte it as function of collision evergy
        let s = normalization * s[j][i];

        let e = newton(
            1e-10,
            s,
            |e| (e + p(e)) / temperature(e) - s,
            |e| e.max(VOID).min(1e4),
        );
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
