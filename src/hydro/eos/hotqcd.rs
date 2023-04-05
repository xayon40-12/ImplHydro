use crate::hydro::HBARC;

const N: usize = 100000;
pub static HOTQCD: [[f64; 4]; N] =
    tofm(unsafe { &std::mem::transmute(*include_bytes!("hrg_hotqcd_eos_binary.dat")) });
// [e, p, s, T]

const fn tofm(arr: &[[f64; 4]; N]) -> [[f64; 4]; N] {
    // let mut arr = *arr;
    // let mut i = 0;
    // while i < N {
    //     arr[i][0] *= HBARC;
    //     i += 1;
    // }
    // arr

    *arr
}

pub fn p(e: f64) -> f64 {
    let e = e * HBARC;
    let e0 = HOTQCD[0][0];
    let e1 = HOTQCD[1][0];
    let e1e0 = e1 - e0;
    let p = if e < e0 {
        let p0 = HOTQCD[0][1];
        p0 * e / e0
    } else {
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let pi = HOTQCD[i][1];
        let pi1 = HOTQCD[i + 1][1];
        pi + (pi1 - pi) * (a - i as f64)
    };
    p / HBARC
}

pub fn dpde(e: f64) -> f64 {
    let eps = 1e-10;
    let dpde = (p(e + eps) - p(e)) / eps;
    dpde
}

pub fn s(e: f64) -> f64 {
    let e = e * HBARC;
    let e0 = HOTQCD[0][0];
    let e1 = HOTQCD[1][0];
    let e1e0 = e1 - e0;
    let s = if e < e0 {
        let s0 = HOTQCD[0][2];
        s0 * e / e0
    } else {
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let si = HOTQCD[i][2];
        let si1 = HOTQCD[i + 1][2];
        si + (si1 - si) * (a - i as f64)
    };
    s
}

#[allow(non_snake_case)]
pub fn T(e: f64) -> f64 {
    let e = e * HBARC;
    let e0 = HOTQCD[0][0];
    let e1 = HOTQCD[1][0];
    let e1e0 = e1 - e0;
    let t = if e < e0 {
        let t0 = HOTQCD[0][3];
        t0 * e / e0
    } else {
        let a = (e - e0) / e1e0;
        let i = (a as usize).min(N - 2);
        let ti = HOTQCD[i][3];
        let ti1 = HOTQCD[i + 1][3];
        ti + (ti1 - ti) * (a - i as f64)
    };
    t / HBARC
}
