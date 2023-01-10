use crate::solver::{
    context::{Boundary, Context},
    run,
    space::kt::{kt, Dir},
    time::{newton::newton, schemes::Scheme},
    utils::ghost,
    Transform,
};

use super::{solve_v, Pressure};

fn constraints(_t: f64, [t00, t01, t02]: [f64; 3]) -> [f64; 3] {
    let m = (t01 * t01 + t02 * t02).sqrt();
    let t00 = t00.max(m);
    [t00, t01, t02]
}

fn gen_transform<'a>(
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
    opt: &Coordinate,
) -> Box<dyn Fn(f64, [f64; 3]) -> [f64; 6] + 'a + Sync> {
    let st = match opt {
        Coordinate::Cartesian => |_| 1.0,
        Coordinate::Milne => |t| t,
    };
    Box::new(move |t, [t00, t01, t02]| {
        let t = st(t);
        let t00 = t00 / t;
        let t01 = t01 / t;
        let t02 = t02 / t;
        let m = (t01 * t01 + t02 * t02).sqrt();
        let sv = solve_v(t00, m, p);
        let v = newton(er, 0.5, |v| sv(v) - v);
        let v = v.max(0.0).min(1.0);
        let e = (t00 - m * v).max(1e-100);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();
        [e, pe, dpde(e), ut, ux, uy]
    })
}

fn eigenvaluesx(_t: f64, [_e, _pe, dpde, ut, ux, _uy]: [f64; 6]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy(t: f64, [e, pe, dpde, ut, ux, uy]: [f64; 6]) -> f64 {
    eigenvaluesx(t, [e, pe, dpde, ut, uy, ux])
}

pub fn f00(t: f64, [e, pe, _, ut, _, _]: [f64; 6]) -> f64 {
    t * ((e + pe) * ut * ut - pe)
}
pub fn f01(t: f64, [e, pe, _, ut, ux, _]: [f64; 6]) -> f64 {
    t * ((e + pe) * ut * ux)
}
pub fn f02(t: f64, [e, pe, _, ut, _, uy]: [f64; 6]) -> f64 {
    t * ((e + pe) * uy * ut)
}

fn f11(t: f64, [e, pe, _, _, ux, _]: [f64; 6]) -> f64 {
    t * ((e + pe) * ux * ux + pe)
}
fn f12(t: f64, [e, pe, _, _, ux, uy]: [f64; 6]) -> f64 {
    t * ((e + pe) * ux * uy)
}

fn f21(t: f64, [e, pe, _, _, ux, uy]: [f64; 6]) -> f64 {
    t * ((e + pe) * uy * ux)
}
fn f22(t: f64, [e, pe, _, _, _, uy]: [f64; 6]) -> f64 {
    t * ((e + pe) * uy * uy + pe)
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 3]; V]; V]; 2],
    constraints: Transform<3, 3>,
    transform: Transform<3, 6>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    opt: &Coordinate,
) -> [f64; 3] {
    let theta = 1.1;

    let t = match opt {
        Coordinate::Cartesian => 1.0,
        Coordinate::Milne => t,
    };

    let pre = &|_t: f64, vs: [f64; 3]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        [m, t01, t02]
    };
    let post = &|_t: f64, vs: [f64; 3]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let t00 = (m * m + k).sqrt();
        [t00, t01, t02]
    };
    let divf1 = kt(
        vs,
        bound,
        pos,
        Dir::X,
        t,
        [&f01, &f11, &f12],
        constraints,
        transform,
        &eigenvaluesx,
        pre,
        post,
        dx,
        theta,
    );
    let divf2 = kt(
        vs,
        bound,
        pos,
        Dir::Y,
        t,
        [&f02, &f21, &f22],
        constraints,
        transform,
        &eigenvaluesy,
        pre,
        post,
        dx,
        theta,
    );

    let s: f64 = match opt {
        Coordinate::Cartesian => 0.0,
        Coordinate::Milne => {
            let [_e, pe, _dpde, _ut, _ux, _uy] =
                transform(t, vs[bound[1](pos[1], V)][bound[0](pos[0], V)]);
            pe
        }
    };
    [
        -divf1[0] - divf2[0] - s,
        -divf1[1] - divf2[1],
        -divf1[2] - divf2[2],
    ]
}
fn flux_exponential<const V: usize>(
    [_ov, vs]: [&[[[f64; 3]; V]; V]; 2],
    _constraints: Transform<3, 3>,
    _transform: Transform<3, 6>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    _dx: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &Coordinate,
) -> [f64; 3] {
    let x = bound[0](pos[0], V);
    let y = bound[1](pos[1], V);
    [-vs[y][x][0], -vs[y][x][1], -vs[y][x][2]]
}

pub type Init2D<'a> = &'a dyn Fn((usize, usize), (f64, f64)) -> [f64; 3];

pub fn hydro2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    opt: Coordinate,
    p: Pressure,
    dpde: Pressure,
    init: Init2D,
    use_exponential: bool,
) -> ([[[f64; 3]; V]; V], f64, usize, usize) {
    let schemename = r.name;
    let mut vs = [[[0.0; 3]; V]; V];
    let names = (["t00", "t01", "t02"], ["e", "pe", "dpde", "ut", "ux", "uy"]);
    let k = [[[[0.0; 3]; V]; V]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init((i, j), (x, y));
        }
    }
    let transform = gen_transform(er, &p, &dpde, &opt);
    let integration = r.integration;
    let flux = if use_exponential {
        flux_exponential
    } else {
        flux
    };
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        transform: &transform,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        local_interaction: [1, 1],   // use a distance of 0 to emulate 1D
        vs,
        k,
        r,
        dt: 1e10,
        dx,
        maxdt,
        er,
        t,
        t0: t,
        tend,
        opt,
        p,
        dpde,
    };
    run(context, name, &schemename, integration, &names)
}
