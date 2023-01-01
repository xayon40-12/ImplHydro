use crate::solver::{
    context::{Boundary, Context, Integration, ToCompute},
    kt::{kt, Dir},
    newton::newton,
    run,
    utils::{ghost, noboundary},
    Constraints,
};

use super::{solve_v, Pressure};

fn gen_constraints<'a>(
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn([f64; 3]) -> [f64; 8] + 'a + Sync> {
    Box::new(move |[t00, t01, v]| {
        let m = t01.abs();
        let t00 = t00.max(m);
        let e = (t00 - m * v).max(1e-100);
        let pe = p(e);
        let v = v.max(0.0).min(1.0);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux).sqrt();
        [t00, t01, v, e, pe, dpde(e), ut, ux]
    })
}

fn eigenvalues([_t00, _t01, _, _e, _pe, dpde, ut, ux]: [f64; 8]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}

pub fn f00([_, _, _, e, pe, _, ut, _]: [f64; 8]) -> f64 {
    (e + pe) * ut * ut - pe
}
pub fn f01([_, _, _, e, pe, _, ut, ux]: [f64; 8]) -> f64 {
    (e + pe) * ut * ux
}
fn f11([_, _, _, e, pe, _, _, ux]: [f64; 8]) -> f64 {
    (e + pe) * ux * ux + pe
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 3]; V]; 1]; 2],
    constraints: Constraints<3, 8>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    er: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
    tocomp: ToCompute,
) -> [f64; 3] {
    let [t00, t01, v, _e, _pe, _dpde, _ut, _ux] = constraints(vs[0][bound[0](pos[0], V)]);
    let rt0 = if tocomp.integrated() || tocomp.all() {
        let theta = 1.5;

        let divf0 = kt(
            vs,
            bound,
            pos,
            Dir::X,
            [&f01, &f11],
            &constraints,
            &eigenvalues,
            dx,
            theta,
        );

        [-divf0[0], -divf0[1]]
    } else {
        [0.0, 0.0]
    };
    let m = t01.abs();
    let sv = solve_v(t00, m);
    let rv = match tocomp {
        ToCompute::All => sv(v),
        ToCompute::NonIntegrated => newton(er, v, |v| sv(v) - v),
        ToCompute::Integrated => 0.0,
    };
    [rt0[0], rt0[1], rv]
}

pub fn hydro1d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: ([[f64; S]; S], Option<[f64; S]>),
    integration: Integration,
    p: Pressure,
    dpde: Pressure,
    init: impl Fn(f64) -> [f64; 3],
) -> ([[[f64; 3]; V]; 1], f64, usize, usize) {
    let mut vs = [[[0.0; 3]; V]];
    let names = ["t00", "t01", "v", "e", "pe", "dpde", "ut", "ux"];
    let mut k = [[[[0.0; 3]; V]]; S];
    let integrated = [true, true, false];
    let v2 = ((V - 1) as f64) / 2.0;
    for i in 0..V {
        let x = (i as f64 - v2) * dx;
        vs[0][i] = init(x);
        for s in 0..S {
            k[s][0][i][2] = vs[0][i][2];
        }
    }
    let constraints = gen_constraints(&p, &dpde);
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &[&ghost, &noboundary], // use noboundary to emulate 1D
        local_interaction: [2, 0],        // use a distance of 0 to emulate 1D
        vs,
        k,
        integrated,
        r,
        dt: 1e10,
        dx,
        maxdt,
        er,
        t,
        tend,
        opt: (),
    };
    run(context, name, integration, &names, &constraints)
}
