use crate::{
    solver::context::Integration::{self, *},
    FLOAT,
};

#[derive(Debug, Clone, Copy)]
pub struct Scheme<const S: usize> {
    pub aij: [[FLOAT; S]; S],
    pub bj: Option<[FLOAT; S]>,
    pub name: &'static str,
    pub integration: Integration,
    pub order: i32,
}

pub fn euler() -> Scheme<1> {
    Scheme {
        aij: [[1.0]],
        bj: None,
        name: "Euler",
        integration: Explicit,
        order: 1,
    }
}
pub fn implicit_euler() -> Scheme<1> {
    Scheme {
        aij: [[1.0]],
        bj: None,
        name: "Euler",
        integration: FixPoint,
        order: 1,
    }
}
pub fn heun() -> Scheme<2> {
    Scheme {
        aij: [[1.0, 0.0], [0.5, 0.5]],
        bj: None,
        name: "Heun",
        integration: Explicit,
        order: 2,
    }
}

pub fn midpoint() -> Scheme<2> {
    Scheme {
        aij: [[0.5, 0.0], [0.0, 1.0]],
        bj: None,
        name: "midpoint",
        integration: Explicit,
        order: 2,
    }
}

pub fn rk4() -> Scheme<4> {
    Scheme {
        aij: [
            [0.5, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0],
        ],
        bj: None,
        name: "rk4",
        integration: Explicit,
        order: 4,
    }
}
pub fn crank_nicolson() -> Scheme<2> {
    Scheme {
        aij: [[0.0, 0.0], [0.5, 0.5]],
        bj: None,
        name: "CrankNicolson",
        integration: FixPoint,
        order: 2,
    }
}
pub fn lobatto_iiic() -> Scheme<2> {
    Scheme {
        aij: [[0.5, -0.5], [0.5, 0.5]],
        bj: None,
        name: "LobattoIIIC",
        integration: FixPoint,
        order: 2,
    }
}
pub fn pareschi() -> Scheme<2> {
    let x = 1.0 - (2.0 as FLOAT).sqrt() / 2.0;
    Scheme {
        aij: [[x, 0.0], [1.0 - 2.0 * x, x]],
        bj: Some([0.5, 0.5]),
        name: "Pareschi",
        integration: FixPoint,
        order: 2,
    }
}
pub fn radauiia2() -> Scheme<2> {
    Scheme {
        aij: [[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]],
        bj: None,
        name: "RadauIIA2",
        integration: FixPoint,
        order: 3,
    }
}
pub fn gauss_legendre_1() -> Scheme<1> {
    Scheme {
        aij: [[0.5]],
        bj: Some([1.0]),
        name: "GL1",
        integration: FixPoint,
        order: 2,
    }
}
pub fn gauss_legendre_2() -> Scheme<2> {
    let sq3 = (3.0 as FLOAT).sqrt();
    Scheme {
        aij: [
            [1.0 / 4.0, 1.0 / 4.0 - 1.0 / 6.0 * sq3],
            [1.0 / 4.0 + 1.0 / 6.0 * sq3, 1.0 / 4.0],
        ],
        bj: Some([1.0 / 2.0, 1.0 / 2.0]),
        name: "GL2",
        integration: FixPoint,
        order: 4,
    }
}

pub fn gauss_legendre_3() -> Scheme<3> {
    let sq15 = (15.0 as FLOAT).sqrt();
    Scheme {
        aij: [
            [
                5.0 / 36.0,
                2.0 / 9.0 - sq15 / 15.0,
                5.0 / 36.0 - sq15 / 30.0,
            ],
            [
                5.0 / 36.0 + sq15 / 24.0,
                2.0 / 9.0,
                5.0 / 36.0 - sq15 / 24.0,
            ],
            [
                5.0 / 36.0 + sq15 / 30.0,
                2.0 / 9.0 + sq15 / 15.0,
                5.0 / 36.0,
            ],
        ],
        bj: Some([5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0]),
        name: "GL3",
        integration: FixPoint,
        order: 6,
    }
}
