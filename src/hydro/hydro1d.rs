use crate::solver::{
    context::{Boundary, Context},
    run,
    space::kt::{kt, Dir},
    time::{newton::newton, schemes::Scheme},
    utils::{ghost, zero},
    Transform,
};

use super::{solve_v, Pressure};

fn constraints(_t: f64, [t00, t01]: [f64; 2]) -> [f64; 2] {
    let m = t01.abs();
    let t00 = t00.max(m);
    [t00, t01]
}
fn gen_transform<'a>(
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn(f64, [f64; 2]) -> [f64; 5] + 'a + Sync> {
    Box::new(move |_t, [t00, t01]| {
        let m = t01.abs();
        let sv = solve_v(t00, m, p);
        let v = newton(er, 0.5, |v| sv(v) - v);
        let v = v.max(0.0).min(1.0);
        let e = (t00 - m * v).max(1e-100);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux).sqrt();
        [e, pe, dpde(e), ut, ux]
    })
}

fn eigenvalues(_t: f64, [_e, _pe, dpde, ut, ux]: [f64; 5]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}

pub fn f00(_t: f64, [e, pe, _, ut, _]: [f64; 5]) -> f64 {
    (e + pe) * ut * ut - pe
}
pub fn f01(_t: f64, [e, pe, _, ut, ux]: [f64; 5]) -> f64 {
    (e + pe) * ut * ux
}
fn f11(_t: f64, [e, pe, _, _, ux]: [f64; 5]) -> f64 {
    (e + pe) * ux * ux + pe
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 2]; V]; 1]; 2],
    constraints: Transform<2, 2>,
    transform: Transform<2, 5>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
) -> [f64; 2] {
    let theta = 1.1;

    let pre = &|_t: f64, vs: [f64; 2]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let k = t01 * t01;
        let m = (t00 * t00 - k).sqrt();
        [m, t01]
    };
    let post = &|_t: f64, vs: [f64; 2]| {
        let m = vs[0];
        let t01 = vs[1];
        let k = t01 * t01;
        let t00 = (m * m + k).sqrt();
        [t00, t01]
    };
    let divf0 = kt(
        vs,
        bound,
        pos,
        Dir::X,
        1.0,
        [&f01, &f11],
        constraints,
        transform,
        &eigenvalues,
        pre,
        post,
        dx,
        theta,
    );

    [-divf0[0], -divf0[1]]
}

pub type Init1D<'a> = &'a dyn Fn(usize, f64) -> [f64; 2];

pub fn hydro1d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Pressure,
    dpde: Pressure,
    init: Init1D,
) -> ([[[f64; 2]; V]; 1], f64, usize, usize) {
    let schemename = r.name;
    let mut vs = [[[0.0; 2]; V]];
    let names = (["t00", "t01"], ["e", "pe", "dpde", "ut", "ux"]);
    let k = [[[[0.0; 2]; V]]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    for i in 0..V {
        let x = (i as f64 - v2) * dx;
        vs[0][i] = init(i, x);
    }
    let transform = gen_transform(er, &p, &dpde);
    let integration = r.integration;
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        transform: &transform,
        boundary: &[&ghost, &zero], // use noboundary to emulate 1D
        local_interaction: [1, 0],  // use a distance of 0 to emulate 1D
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
        opt: (),
        p,
        dpde,
    };
    run(context, name, &schemename, integration, &names)
}
