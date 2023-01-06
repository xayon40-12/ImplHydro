use crate::hydro::hydro1d::{f00, f01};

use super::Pressure;

pub fn init_riemann<'a>(
    p: Pressure<'a>,
    dpde: Pressure<'a>,
    use_void: bool,
) -> Box<dyn Fn(f64) -> [f64; 2] + 'a> {
    let el = 10.0;
    let er = if use_void { 1e-100 } else { 1.0 };
    Box::new(move |x| {
        // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
        let e = if x < 0.0 { el } else { er };
        let vars = [e, p(e), dpde(e), 1.0, 0.0];
        [f00(vars), f01(vars)]
    })
}
