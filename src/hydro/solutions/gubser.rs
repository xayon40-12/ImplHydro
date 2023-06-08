use crate::hydro::ideal::ideal2d::f0;

use crate::hydro::Eos;

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

pub fn init_gubser<'a>(
    t0: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; 3] + 'a> {
    Box::new(move |_, (x, y)| {
        let [e, ut, ux, uy] = gubser(x, y, t0);
        let vars = [e, p(e), dpde(e), ut, ux, uy];
        f0(t0, vars)
    })
}
