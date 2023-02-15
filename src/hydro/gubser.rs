use crate::hydro::hydro2d::{f00, f01, f02};

use super::Eos;

pub fn gubser(x: f64, y: f64, t: f64) -> [f64; 4] {
    let r = (x * x + y * y).sqrt();
    let e = 2.0f64.powf(8.0 / 3.0)
        / (t.powf(4.0 / 3.0)
            * (1.0 + 2.0 * (t * t + r * r) + (t * t - r * r).powf(2.0)).powf(4.0 / 3.0));
    let k = (2.0 * t * r / (1.0 + t * t + r * r)).atanh();
    let ut = k.cosh();
    let ux;
    let uy;
    if r == 0.0 {
        ux = 0.0;
        uy = 0.0;
    } else {
        ux = x / r * k.sinh();
        uy = y / r * k.sinh();
    }
    [e, ut, ux, uy]
}
pub fn gubser_err<const V: usize>(v: [[[f64; 4]; V]; V], t: f64, dx: f64, p: Eos) -> [f64; 2] {
    let v2 = ((V - 1) as f64) / 2.0;
    let mut maxerrt00 = 0.0f64;
    let mut meanerrt00 = 0.0f64;
    for i in 0..V {
        for j in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            let [e, ut, _, _] = gubser(x, y, t);

            let vt00 = v[j][i][0];
            let gt00 = (e + p(e)) * ut * ut - p(e);
            let errt00 = (vt00 - gt00).abs() / vt00.abs().max(gt00.abs());
            maxerrt00 = maxerrt00.max(errt00);
            meanerrt00 += errt00;
        }
    }
    meanerrt00 /= (V * V) as f64;
    [maxerrt00, meanerrt00]
}

pub fn init_gubser<'a>(
    t0: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; 3] + 'a> {
    Box::new(move |_, (x, y)| {
        let [e, ut, ux, uy] = gubser(x, y, t0);
        let vars = [e, p(e), dpde(e), ut, ux, uy];
        [f00(t0, vars), f01(t0, vars), f02(t0, vars)]
    })
}
