use crate::hydro::{eos::cubic_spline::cubic_spline_0, HBARC};
use lazy_static::lazy_static;

use super::cubic_spline::{cubic, cubic_diff};

const N: usize = 100000;

// [e, p, s, T]
// [GeV.fm^-3, GeV.fm^-3, fm^-3, GeV]
static HOTQCD_GEV: &[[f64; 4]; N] =
    unsafe { &std::mem::transmute(*include_bytes!("hrg_hotqcd_eos_binary.dat")) };

static HOTQCD_ID_E: usize = 0;
static HOTQCD_ID_P: usize = 1;
// static HOTQCD_ID_S: usize = 2;
static HOTQCD_ID_T: usize = 3;

lazy_static! {
    // [e, p, s, T]
    // [fm^-4, fm^-4, fm^-3, fm^-1]
    static ref HOTQCD_FM: Box<[[f64; 4]; N]> = {
        let mut res: Box<[[f64; 4]; N]> = Box::new(*HOTQCD_GEV);
        // before [GeV.fm^-3, GeV.fm^-3, fm^-3, GeV]
        for i in 0..N {
            res[i][0] /= HBARC;
            res[i][1] /= HBARC;
            res[i][3] /= HBARC;
        }
        // after [fm^-4, fm^-4, fm^-3, fm^-1]
        res
    };
    static ref HOTQCD_P: Box<[[f64; 4]; N]> = cubic_spline_0(HOTQCD_ID_E, HOTQCD_ID_P, &HOTQCD_FM);
    static ref HOTQCD_T: Box<[[f64; 4]; N]> = cubic_spline_0(HOTQCD_ID_E, HOTQCD_ID_T, &HOTQCD_FM);
}

fn hot_cubic(e: f64, spline: &[[f64; 4]; N], cubic: &dyn Fn(f64, &[f64; 4], f64) -> f64) -> f64 {
    let e0 = HOTQCD_FM[0][0];
    let e1 = HOTQCD_FM[1][0];
    let e1e0 = e1 - e0;
    let h = if e < e0 {
        let a = e0 / e1e0;
        cubic(e / e0, &spline[0], a)
    } else {
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let x = a - i as f64;
        cubic(x, &spline[i + 1], 1.0)
    };
    h
}

pub fn p(e: f64) -> f64 {
    hot_cubic(e, &HOTQCD_P, &cubic)
}

pub fn dpde(e: f64) -> f64 {
    hot_cubic(e, &HOTQCD_P, &cubic_diff)
}

#[allow(non_snake_case)]
pub fn T(e: f64) -> f64 {
    hot_cubic(e, &HOTQCD_T, &cubic)
}

#[test]
pub fn test_hotqcd() {
    let n = 10000;
    for i in 0..n {
        let x = i as f64 / n as f64;
        let e = x * 1e0;
        let p = p(e);
        let dpde = dpde(e);
        println!("@ {e:.5e} {p:.5e} {dpde:.5e}");
    }
}
