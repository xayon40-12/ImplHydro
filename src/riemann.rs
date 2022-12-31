use crate::hydro1d::{f00, f01};

pub fn init_riemann() -> Box<dyn Fn(f64) -> [f64; 3]> {
    Box::new(move |x| {
        // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
        let e = if x < 0.0 { 10.0 } else { 1e-100 };
        let vars = [0.0, 0.0, e, 1.0, 0.0];
        [f00(vars), f01(vars), e]
    })
}
