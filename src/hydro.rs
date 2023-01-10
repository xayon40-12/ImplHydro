pub mod eos;
pub mod gubser;
pub mod hydro1d;
pub mod hydro2d;
pub mod riemann;

pub type Pressure<'a> = &'a (dyn Fn(f64) -> f64 + Sync);

pub mod ideal_gas {
    pub fn p(e: f64) -> f64 {
        e / 3.0
    }

    pub fn dpde(_e: f64) -> f64 {
        1.0 / 3.0
    }
}

pub fn solve_v<'a>(t00: f64, m: f64, p: Pressure<'a>) -> Box<dyn Fn(f64) -> f64 + 'a> {
    Box::new(move |v| {
        let e = (t00 - m * v).max(1e-100);
        let v = m / (t00 + p(e));
        v
    })
}
