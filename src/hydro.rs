use crate::solver::context::BArr;

pub mod eos;
pub mod ideal;
pub mod isosurface;
pub mod solutions;
pub mod utils;
pub mod viscous;

pub type Eos<'a> = &'a (dyn Fn(f64) -> f64 + Sync);

pub type Init1D<'a, const F: usize> = &'a dyn Fn(usize, f64) -> [f64; F];
pub type Init2D<'a, const F: usize> = &'a dyn Fn((usize, usize), (f64, f64)) -> [f64; F];
pub type Init3D<'a, const F: usize> =
    &'a dyn Fn((usize, usize, usize), (f64, f64, f64)) -> [f64; F];

pub const VOID: f64 = 1e-100;
pub const HBARC: f64 = 0.1973; // GeV.fm

#[derive(Debug, Clone, Copy)]
pub enum Dim {
    D1,
    D2,
    D3,
}

impl Dim {
    pub const fn value(&self) -> usize {
        match self {
            Dim::D1 => 1,
            Dim::D2 => 2,
            Dim::D3 => 3,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Viscosity {
    Ideal,
    Bulk(f64, f64),                              // (bulk,energycut)
    Shear((f64, f64, f64), f64),                 // ((eta/s_min, eta/s_slope, eta/s_crv),energycut)
    Both((f64, f64, f64), (f64, f64, f64), f64), // ((eta/s_min, eta/s_slope, eta/s_crv),(zeta/s_max,zeta/s_width,zeta/s_peak),energycut)
}

impl Viscosity {
    pub const fn nb_fields(&self, dim: usize) -> usize {
        let d = dim + 1;
        let bulk = 1;
        let shear = d * (d - 1) / 2;
        match self {
            Viscosity::Ideal => d,
            Viscosity::Bulk(_, _) => d + bulk,
            Viscosity::Shear(_, _) => d + shear,
            Viscosity::Both(_, _, _) => d + shear + bulk,
        }
    }
    pub const fn nb_transforms(&self, dim: usize) -> usize {
        let f = self.nb_fields(dim) + 3;
        let d = dim + 1;
        match self {
            Viscosity::Ideal => f,
            Viscosity::Bulk(_, _) => f,
            Viscosity::Shear(_, _) => f + d,   // +d for pi00..pi02
            Viscosity::Both(_, _, _) => f + d, // +d for pi00..pi02
        }
    }
}

pub type HydroOutput<
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const F: usize,
    const C: usize,
> = Option<(
    (BArr<F, VX, VY, VZ>, BArr<C, VX, VY, VZ>),
    f64,
    usize,
    usize,
)>;

pub const FREESTREAM_2D: usize = 11;

pub const F_IDEAL_1D: usize = Viscosity::Ideal.nb_fields(Dim::D1.value());
pub const F_BULK_1D: usize = Viscosity::Bulk(0.0, 0.0).nb_fields(Dim::D1.value());
pub const F_SHEAR_1D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D1.value());
pub const F_BOTH_1D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D1.value());
pub const F_IDEAL_2D: usize = Viscosity::Ideal.nb_fields(Dim::D2.value());
pub const F_BULK_2D: usize = Viscosity::Bulk(0.0, 0.0).nb_fields(Dim::D2.value());
pub const F_SHEAR_2D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D2.value());
pub const F_BOTH_2D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D2.value());
pub const F_IDEAL_3D: usize = Viscosity::Ideal.nb_fields(Dim::D3.value());
pub const F_BULK_3D: usize = Viscosity::Bulk(0.0, 0.0).nb_fields(Dim::D3.value());
pub const F_SHEAR_3D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D3.value()) - 1;
pub const F_BOTH_3D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_fields(Dim::D3.value()) - 1;

pub const C_IDEAL_1D: usize = Viscosity::Ideal.nb_transforms(Dim::D1.value());
pub const C_BULK_1D: usize = Viscosity::Bulk(0.0, 0.0).nb_transforms(Dim::D1.value());
pub const C_SHEAR_1D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D1.value());
pub const C_BOTH_1D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D1.value());
pub const C_IDEAL_2D: usize = Viscosity::Ideal.nb_transforms(Dim::D2.value());
pub const C_BULK_2D: usize = Viscosity::Bulk(0.0, 0.0).nb_transforms(Dim::D2.value());
pub const C_SHEAR_2D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D2.value());
pub const C_BOTH_2D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D2.value());
pub const C_MILNE_SHEAR_2D: usize =
    Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D2.value()) + 1; // +1 for pi33
pub const C_MILNE_BOTH_2D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D2.value()) + 1; // +1 for pi33
pub const C_IDEAL_3D: usize = Viscosity::Ideal.nb_transforms(Dim::D3.value());
pub const C_BULK_3D: usize = Viscosity::Bulk(0.0, 0.0).nb_transforms(Dim::D3.value());
pub const C_SHEAR_3D: usize = Viscosity::Shear((0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D3.value());
pub const C_BOTH_3D: usize =
    Viscosity::Both((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0).nb_transforms(Dim::D3.value());

pub mod ideal_gas {
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
}

pub fn solve_v<'a>(t00: f64, m: f64, p: Eos<'a>) -> Box<dyn Fn(f64) -> f64 + 'a> {
    Box::new(move |v| {
        let e = (t00 - m * v).max(VOID);
        let v = m / (t00 + p(e));
        v
    })
}
