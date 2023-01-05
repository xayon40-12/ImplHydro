use super::context::Integration::{self, *};

#[derive(Debug, Clone, Copy)]
pub struct Scheme<const S: usize> {
    pub aij: [[f64; S]; S],
    pub bj: Option<[f64; S]>,
    pub name: &'static str,
    pub integration: Integration,
}

pub fn euler() -> Scheme<1> {
    Scheme {
        aij: [[1.0]],
        bj: None,
        name: "Euler",
        integration: Explicit,
    }
}
pub fn implicit_euler() -> Scheme<1> {
    Scheme {
        aij: [[1.0]],
        bj: None,
        name: "Euler",
        integration: FixPoint,
    }
}
pub fn heun() -> Scheme<2> {
    Scheme {
        aij: [[1.0, 0.0], [0.5, 0.5]],
        bj: None,
        name: "Heun",
        integration: Explicit,
    }
}
pub fn radauiia2() -> Scheme<2> {
    Scheme {
        aij: [[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]],
        bj: None,
        name: "RadauIIA2",
        integration: FixPoint,
    }
}
pub fn gauss_legendre_1() -> Scheme<1> {
    Scheme {
        aij: [[0.5]],
        bj: Some([1.0]),
        name: "GL1",
        integration: FixPoint,
    }
}
pub fn gauss_legendre_2() -> Scheme<2> {
    let sq3 = 3.0f64.sqrt();
    Scheme {
        aij: [
            [1.0 / 4.0, 1.0 / 4.0 - 1.0 / 6.0 * sq3],
            [1.0 / 4.0 + 1.0 / 6.0 * sq3, 1.0 / 4.0],
        ],
        bj: Some([1.0 / 2.0, 1.0 / 2.0]),
        name: "GL2",
        integration: FixPoint,
    }
}
