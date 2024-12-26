pub mod cubic_spline;

pub mod wb;

pub enum EOSs {
    ConformalMassless,
    WB,
}

pub mod conformal_massless {
    use std::f64::consts::PI;

    use crate::FLOAT;

    pub fn p(e: FLOAT) -> FLOAT {
        e / 3.0
    }

    pub fn dpde(_e: FLOAT) -> FLOAT {
        1.0 / 3.0
    }

    #[allow(non_snake_case)]
    pub fn T(e: FLOAT) -> FLOAT {
        let pi = PI as FLOAT;
        (30.0 * e / (pi * pi * (16.0 + 21.0 / 2.0))).powf(0.25)
    }

    pub fn s(e: FLOAT) -> FLOAT {
        (e + p(e)) / (T(e) + 1e-100)
    }
}
