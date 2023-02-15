use crate::hydro::ideal::ideal1d::{f00, f01};

use crate::hydro::{Eos, VOID};

pub fn init_riemann<'a>(
    t0: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
    use_void: bool,
) -> Box<dyn Fn(usize, f64) -> [f64; 2] + 'a> {
    let el = 10.0;
    let er = if use_void { VOID } else { 1.0 };
    Box::new(move |_, x| {
        let e = if x < 0.0 { el } else { er };
        let vars = [e, p(e), dpde(e), 1.0, 0.0];
        [f00(t0, vars), f01(t0, vars)]
    })
}
