use crate::hydro::ideal::ideal1d::f0;

use crate::hydro::{Eos, VOID};
use crate::FLOAT;

pub fn init_riemann<'a>(
    t0: FLOAT,
    p: Eos<'a>,
    dpde: Eos<'a>,
    use_void: bool,
) -> Box<dyn Fn(usize, FLOAT) -> [FLOAT; 2] + 'a> {
    let el = 10.0;
    let er = if use_void { VOID } else { 1.0 };
    Box::new(move |_, x| {
        let e = if x < 0.0 { el } else { er };
        let vars = [e, p(e), dpde(e), 1.0, 0.0];
        f0(t0, vars)
    })
}
