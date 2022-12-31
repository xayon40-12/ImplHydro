use crate::hydro2d::{f00, f01, f02, p};

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
pub fn gubser_err<const V: usize>(v: [[[f64; 4]; V]; V], t: f64, dx: f64) -> [f64; 4] {
    let v2 = ((V - 1) as f64) / 2.0;
    let mut maxerr = 0.0f64;
    let mut meanerr = 0.0f64;
    let mut maxerrt00 = 0.0f64;
    let mut meanerrt00 = 0.0f64;
    for i in 0..V {
        for j in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            let [e, ut, _, _] = gubser(x, y, t);

            let ve = v[j][i][3];
            let err = (ve - e).abs() / ve.abs().max(e.abs());
            maxerr = maxerr.max(err);
            meanerr += err;

            let vt00 = v[j][i][0];
            let gt00 = (e + p(e)) * ut * ut - p(e);
            let errt00 = (vt00 - gt00).abs() / vt00.abs().max(gt00.abs());
            maxerrt00 = maxerrt00.max(errt00);
            meanerrt00 += errt00;
        }
    }
    meanerr /= (V * V) as f64;
    meanerrt00 /= (V * V) as f64;
    [maxerr, meanerr, maxerrt00, meanerrt00]
}

pub fn init_gubser(t0: f64) -> Box<dyn Fn(f64, f64) -> [f64; 4]> {
    Box::new(move |x, y| {
        // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
        let [e, ut, ux, uy] = gubser(x, y, t0);
        let vars = [0.0, 0.0, 0.0, e, ut, ux, uy];
        [f00(vars), f01(vars), f02(vars), e]
    })
}
