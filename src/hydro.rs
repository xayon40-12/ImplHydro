pub mod eos;
pub mod from_file;
pub mod gubser;
pub mod hydro1d;
pub mod hydro2d;
pub mod riemann;
pub mod viscoushydro2d;

pub type Eos<'a> = &'a (dyn Fn(f64) -> f64 + Sync);

pub type Init1D<'a, const F: usize> = &'a dyn Fn(usize, f64) -> [f64; F];
pub type Init2D<'a, const F: usize> = &'a dyn Fn((usize, usize), (f64, f64)) -> [f64; F];

pub static VOID: f64 = 1e-100;

#[derive(Debug, Clone, Copy)]
pub enum Dim {
    D1,
    D2,
}

impl Dim {
    pub const fn value(&self) -> usize {
        match self {
            Dim::D1 => 1,
            Dim::D2 => 2,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Viscosity {
    Ideal,
    Bulk,
    Shear,
    Both,
}

impl Viscosity {
    pub const fn nb_fields(&self, dim: usize) -> usize {
        let d = dim + 1;
        let bulk = 1;
        let shear = d * (d + 1) / 2;
        match self {
            Viscosity::Ideal => d,
            Viscosity::Bulk => d + bulk,
            Viscosity::Shear => d + shear,
            Viscosity::Both => d + shear + bulk,
        }
    }
}

pub type HydroOutput<const VX: usize, const VY: usize, const F: usize> =
    Option<([[[f64; F]; VX]; VY], f64, usize, usize)>;

pub const F_IDEAL_1D: usize = Viscosity::Ideal.nb_fields(Dim::D1.value());
pub const F_BULK_1D: usize = Viscosity::Bulk.nb_fields(Dim::D1.value());
pub const F_SHEAR_1D: usize = Viscosity::Shear.nb_fields(Dim::D1.value());
pub const F_BOTH_1D: usize = Viscosity::Both.nb_fields(Dim::D1.value());
pub const F_IDEAL_2D: usize = Viscosity::Ideal.nb_fields(Dim::D2.value());
pub const F_BULK_2D: usize = Viscosity::Bulk.nb_fields(Dim::D2.value());
pub const F_SHEAR_2D: usize = Viscosity::Shear.nb_fields(Dim::D2.value());
pub const F_BOTH_2D: usize = Viscosity::Both.nb_fields(Dim::D2.value());

pub mod ideal_gas {
    pub fn p(e: f64) -> f64 {
        e / 3.0
    }

    pub fn dpde(_e: f64) -> f64 {
        1.0 / 3.0
    }
}

pub fn solve_v<'a>(t00: f64, m: f64, p: Eos<'a>) -> Box<dyn Fn(f64) -> f64 + 'a> {
    Box::new(move |v| {
        let e = (t00 - m * v).max(VOID);
        let v = m / (t00 + p(e));
        v
    })
}
