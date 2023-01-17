pub mod eos;
pub mod from_file;
pub mod gubser;
pub mod hydro1d;
pub mod hydro2d;
pub mod riemann;
pub mod viscoushydro2d;

pub type Pressure<'a> = &'a (dyn Fn(f64) -> f64 + Sync);

pub type Init2D<'a, const F: usize> = &'a dyn Fn((usize, usize), (f64, f64)) -> [f64; F];

pub static VOID: f64 = 1e-100;

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
        let e = (t00 - m * v).max(VOID);
        let v = m / (t00 + p(e));
        v
    })
}
