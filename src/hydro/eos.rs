pub mod cubic_spline;
pub mod hotqcd;
pub mod wb;

pub enum EOSs {
    ConformalMassless,
    WB,
    HotQCD,
    HotQCDLog,
}

pub mod conformal_massless {
    use std::f64::consts::PI;

    pub fn p(e: f64) -> f64 {
        e / 3.0
    }

    pub fn dpde(_e: f64) -> f64 {
        1.0 / 3.0
    }

    #[allow(non_snake_case)]
    pub fn T(e: f64) -> f64 {
        (30.0 * e / (PI * PI * (16.0 + 21.0 / 2.0))).powf(0.25)
    }

    pub fn s(e: f64) -> f64 {
        (e + p(e)) / (T(e) + 1e-100)
    }
}
