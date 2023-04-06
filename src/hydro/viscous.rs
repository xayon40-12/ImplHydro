use crate::{
    hydro::{Eos, VOID},
    solver::time::newton::newton,
};

use super::{FREESTREAM_2D, F_BOTH_2D};

pub mod viscous2d;

pub fn init_from_entropy_density_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    s: [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
    temperature: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_BOTH_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let s = s[j][i];

        let e = newton(
            1e-10,
            s,
            |e| (e + p(e)) / temperature(e) - s,
            |e| e.max(VOID).min(1e10),
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
            0.0,
        ];
        viscous2d::fitutpi(t0, vars)
    })
}

pub fn init_from_freestream_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    trs: [[[f64; FREESTREAM_2D]; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_BOTH_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let e = trs[j][i][0].max(VOID);
        let ut = trs[j][i][1];
        let ux = trs[j][i][2];
        let uy = trs[j][i][3];
        // TODO use viscosity from freestream
        // let pi00 = trs[j][i][4];
        // let pi01 = trs[j][i][5];
        // let pi02 = trs[j][i][6];
        // let pi11 = trs[j][i][7];
        // let pi12 = trs[j][i][8];
        // let pi22 = trs[j][i][9];
        // let pi33 = (pi00 - pi11 - pi22) / (t0 * t0);
        // let bulk = trs[j][i][10];
        let vars = [
            e,
            p(e),
            dpde(e),
            ut,
            ux,
            uy,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            // pi00,
            // pi01,
            // pi02,
            // pi11,
            // pi12,
            // pi22,
            // pi33,
            // bulk,
        ];
        // println!("vars: {:?}", vars);
        viscous2d::fitutpi(t0, vars)
    })
}
