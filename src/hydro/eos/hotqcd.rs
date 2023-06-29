use crate::{boxarray, hydro::HBARC};
use lazy_static::lazy_static;

const N: usize = 100000;
const N1: usize = 100001;

lazy_static! {
    static ref HOTQCD: Box<[[f64; 4]; N]> = {
        let mut tmp: [[f64; 4]; N] =
            *unsafe { &std::mem::transmute(*include_bytes!("hrg_hotqcd_eos_binary.dat")) };
        // [e, p, s, T]
        // before [GeV.fm^-3, GeV.fm^-3, fm^-3, GeV]
        for i in 0..N {
            tmp[i][0] /= HBARC;
            tmp[i][1] /= HBARC;
            tmp[i][3] /= HBARC;
        }
        // after [fm^-4, fm^-4, fm^-3, fm^-1]
        Box::new(tmp)
    };
    static ref HOTQCD_P: Box<[[f64; 4]; N]> = cubic_spline_0(1, &HOTQCD);
    static ref HOTQCD_T: Box<[[f64; 4]; N]> = cubic_spline_0(3, &HOTQCD);
}

fn cubic_spline_0(id: usize, arr: &[[f64; 4]; N]) -> Box<[[f64; 4]; N]> {
    let mut tmp: Box<[[f64; 3]; N1]> = boxarray(0.0f64);
    // WARNING: as the first segment down to 0 is not the same lenght as the other segments, the spline does not work properly for this first segment
    tmp[0][0] = 0.0;
    tmp[0][2] = 2.0;
    for i in 0..N {
        tmp[i + 1][0] = arr[i][id];
        tmp[i + 1][2] = 4.0;
    }
    tmp[N][2] = 2.0;

    tmp[0][1] = 3.0 * (tmp[1][0] - tmp[0][0]);
    tmp[N][1] = 3.0 * (tmp[N][0] - tmp[N - 1][0]);
    for i in 1..N {
        tmp[i][1] = 3.0 * (tmp[i + 1][0] - tmp[i - 1][0]);
    }

    for i in 1..=N {
        let d = 1.0 / tmp[i - 1][2];
        tmp[i][2] -= d;
        tmp[i][1] -= d * tmp[i - 1][1];
    }
    for i in (0..N).rev() {
        tmp[i + 1][1] /= tmp[i + 1][2];
        tmp[i][1] -= tmp[i + 1][1];
    }
    tmp[0][1] /= tmp[0][2];

    let mut pols: Box<[[f64; 4]; N]> = boxarray(0.0f64);
    for i in 0..N {
        pols[i][0] = tmp[i][0];
        pols[i][1] = tmp[i][1];
        pols[i][2] = 3.0 * (tmp[i + 1][0] - tmp[i][0]) - 2.0 * tmp[i][1] - tmp[i + 1][1];
        pols[i][3] = 2.0 * (tmp[i][0] - tmp[i + 1][0]) + tmp[i][1] + tmp[i + 1][1];
    }

    pols
}

fn cubic(x: f64, p: &[f64; 4]) -> f64 {
    p[0] + x * (p[1] + x * (p[2] + x * p[3]))
}
fn cubic_diff(x: f64, p: &[f64; 4]) -> f64 {
    p[1] + x * (2.0 * p[2] + x * 3.0 * p[3])
}

fn hot_cubic(e: f64, spline: &[[f64; 4]; N], cubic: &dyn Fn(f64, &[f64; 4]) -> f64) -> f64 {
    let e0 = HOTQCD[0][0];
    let h = if e < e0 {
        cubic(e / e0, &spline[0])
    } else {
        let e1 = HOTQCD[1][0];
        let e1e0 = e1 - e0;
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let x = a - i as f64;
        cubic(x, &spline[i + 1])
    };
    h
}

// fn hot(e: f64, c: usize) -> f64 {
//     let e0 = HOTQCD[0][0];
//     let h = if e < e0 {
//         let h0 = HOTQCD[0][c];
//         h0 * e / e0
//     } else {
//         let e1 = HOTQCD[1][0];
//         let e1e0 = e1 - e0;
//         let a = (e - e0) / e1e0;
//         let i = (a as usize).min(N - 2);
//         let hi = HOTQCD[i][c];
//         let hi1 = HOTQCD[i + 1][c];
//         hi + (hi1 - hi) * (a - i as f64)
//     };
//     h
// }

// fn dhot(e: f64, c: usize) -> f64 {
//     let e0 = HOTQCD[0][0];
//     let e1 = HOTQCD[1][0];
//     let e1e0 = e1 - e0;
//     let h = if e < e0 * 0.5 {
//         let h0 = HOTQCD[0][c];
//         let dh = h0 / e0;
//         dh
//     } else if e < (e0 + e1) * 0.5 {
//         let a = (e - e0 * 0.5) / ((e0 + e1) * 0.5 - e0 * 0.5);
//         let h0 = HOTQCD[0][c];
//         let h1 = HOTQCD[1][c];
//         let dh0 = h0 / e0;
//         let dh1 = (h1 - h0) / e1e0;
//         dh0 + (dh1 - dh0) * a
//     } else {
//         let a = (e - (e0 + e1) * 0.5) / e1e0;
//         let i = (a as usize).min(N - 3);
//         let hi = HOTQCD[i][c];
//         let hi1 = HOTQCD[i + 1][c];
//         let hi2 = HOTQCD[i + 2][c];
//         let dhi1 = (hi1 - hi) / e1e0;
//         let dhi2 = (hi2 - hi1) / e1e0;
//         dhi1 + (dhi2 - dhi1) * (a - i as f64)
//     };
//     h
// }

// #[test]
// pub fn test_dhot() {
//     let n = 10000;
//     let m = 1;
//     for i in 0..m * n {
//         let e = i as f64 / n as f64;
//         let t = T(e);
//         let h = hot(e, 1);
//         let d = dhot(e, 1);
//         println!("@ {:e} {:e} {:e} {:e}", e, t, h, d);
//     }
// }

pub fn p(e: f64) -> f64 {
    // hot(e, 1)
    hot_cubic(e, &HOTQCD_P, &cubic)
}

pub fn dpde(e: f64) -> f64 {
    // dhot(e, 1)
    hot_cubic(e, &HOTQCD_P, &cubic_diff)
}

#[allow(non_snake_case)]
pub fn T(e: f64) -> f64 {
    // hot(e, 3)
    hot_cubic(e, &HOTQCD_T, &cubic)
}
