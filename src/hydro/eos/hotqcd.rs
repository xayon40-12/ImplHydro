use crate::hydro::HBARC;

const N: usize = 100000;
pub static HOTQCD: [[f64; 4]; N] =
    *unsafe { &std::mem::transmute(*include_bytes!("hrg_hotqcd_eos_binary.dat")) };
// [e, p, s, T]

fn hot(e: f64, c: usize) -> f64 {
    let e = e * HBARC;
    let e0 = HOTQCD[0][0];
    let e1 = HOTQCD[1][0];
    let e1e0 = e1 - e0;
    let h = if e < e0 {
        let h0 = HOTQCD[0][c];
        h0 * e / e0
    } else {
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let hi = HOTQCD[i][c];
        let hi1 = HOTQCD[i + 1][c];
        hi + (hi1 - hi) * (a - i as f64)
    };
    h / HBARC
}

fn dhot(e: f64, c: usize) -> f64 {
    let e = e * HBARC;
    let e0 = HOTQCD[0][0];
    let e1 = HOTQCD[1][0];
    let e1e0 = e1 - e0;
    let h = if e < e0 * 0.5 {
        let h0 = HOTQCD[0][c];
        let dh = h0 / e0;
        dh
    } else if e < (e0 + e1) * 0.5 {
        let a = (e - e0 * 0.5) / ((e0 + e1) * 0.5 - e0 * 0.5);
        let h0 = HOTQCD[0][c];
        let h1 = HOTQCD[1][c];
        let dh0 = h0 / e0;
        let dh1 = (h1 - h0) / e1e0;
        dh0 + (dh1 - dh0) * a
    } else {
        let a = (e - (e0 + e1) * 0.5) / e1e0;
        let i = (a as usize).min(N - 3);
        let hi = HOTQCD[i][c];
        let hi1 = HOTQCD[i + 1][c];
        let hi2 = HOTQCD[i + 2][c];
        let dhi1 = (hi1 - hi) / e1e0;
        let dhi2 = (hi2 - hi1) / e1e0;
        dhi1 + (dhi2 - dhi1) * (a - i as f64)
    };
    h
}

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
    hot(e, 1)
}

pub fn dpde(e: f64) -> f64 {
    dhot(e, 1)
    // let eps = 1e-10;
    // let dpde = (p(e + eps) - p(e)) / eps;
    // dpde
}

pub fn s(e: f64) -> f64 {
    hot(e, 2)
}

#[allow(non_snake_case)]
pub fn T(e: f64) -> f64 {
    hot(e, 3)
}
