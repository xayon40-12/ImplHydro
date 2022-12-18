use crate::{
    kt::{kt, Dir},
    newton::{Boundary, Context},
    solver::run,
    utils::ghost,
};

fn p(e: f64) -> f64 {
    e / 3.0
}

fn dpde(_e: f64) -> f64 {
    1.0 / 3.0
}

fn constraints([t00, t01, t02, e]: [f64; 4]) -> [f64; 7] {
    let t00 = t00.max((t01 * t01 + t02 * t02).sqrt());
    let e = e.max(1e-100);
    let ut = ((t00 + p(e)) / (e + p(e))).sqrt().max(1.0);
    let ux = t01 / ((e + p(e)) * ut);
    let uy = t02 / ((e + p(e)) * ut);
    let ut = (1.0 + ux * ux + uy * uy).sqrt();
    [t00, t01, t02, e, ut, ux, uy]
}

fn eigenvaluesx([_t00, _t01, _t02, e, ut, ux, _uy]: [f64; 7]) -> f64 {
    let vs2 = dpde(e);
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy([t00, t01, t02, e, ut, ux, uy]: [f64; 7]) -> f64 {
    eigenvaluesx([t00, t01, t02, e, ut, uy, ux])
}

pub fn f00([_, _, _, e, ut, _, _]: [f64; 7]) -> f64 {
    (e + p(e)) * ut * ut - p(e)
}
pub fn f01([_, _, _, e, ut, ux, _]: [f64; 7]) -> f64 {
    (e + p(e)) * ut * ux
}
pub fn f02([_, _, _, e, ut, _, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * uy * ut
}

fn f11([_, _, _, e, _, ux, _]: [f64; 7]) -> f64 {
    (e + p(e)) * ux * ux + p(e)
}
fn f12([_, _, _, e, _, ux, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * ux * uy
}

fn f21([_, _, _, e, _, ux, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * uy * ux
}
fn f22([_, _, _, e, _, _, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * uy * uy + p(e)
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, v]: [&[[[f64; 4]; V]; V]; 2],
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    opt: &Coordinate,
) -> [f64; 4] {
    let theta = 1.5;

    let divf1 = kt(
        v,
        bound,
        pos,
        Dir::X,
        [&f01, &f11, &f12],
        &constraints,
        &eigenvaluesx,
        dx,
        theta,
    );
    let divf2 = kt(
        v,
        bound,
        pos,
        Dir::Y,
        [&f02, &f21, &f22],
        &constraints,
        &eigenvaluesy,
        dx,
        theta,
    );

    let [t00, t01, t02, e, ut, ux, uy] = constraints(v[bound[1](pos[1], V)][bound[0](pos[0], V)]);
    let re = t00 - (t01 * t01 + t02 * t02) / (t00 + p(e));
    let source: [f64; 3] = match opt {
        Coordinate::Cartesian => [0.0; 3],
        Coordinate::Milne => [
            (e + p(e)) * ut * ut / t,
            (e + p(e)) * ut * ux / t,
            (e + p(e)) * ut * uy / t,
        ],
    };
    let rt0 = [
        -divf1[0] - divf2[0] - source[0],
        -divf1[1] - divf2[1] - source[1],
        -divf1[2] - divf2[2] - source[2],
    ];
    [rt0[0], rt0[1], rt0[2], re]
}

pub fn hydro2d<const V: usize>(
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    opt: Coordinate,
    init: impl Fn(f64, f64) -> [f64; 4],
) -> ([[[f64; 4]; V]; V], f64) {
    let mut vs = [[[0.0; 4]; V]; V];
    let names = ["t00", "t01", "t02", "e", "ut", "ux", "uy"];
    let k = [[[[0.0; 4]; V]; V]];
    let integrated = [true, true, true, false];
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init(x, y);
        }
    }
    let context = Context {
        fun: &flux,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        local_interaction: [2, 2],   // use a distance of 0 to emulate 1D
        vs,
        k,
        integrated,
        r: [[1.0]],
        dt: 1e10,
        dx,
        maxdt,
        er,
        t,
        tend,
        opt,
    };
    let (vals, t, cost, tsteps) = run(context, &names, &constraints);
    println!("cost: {}, tsteps: {}", cost, tsteps);
    (vals, t)
}
