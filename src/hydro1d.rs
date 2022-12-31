use crate::{
    context::{Boundary, Context, Integration, ToCompute},
    kt::{kt, Dir},
    newton::newton,
    solver::run,
    utils::{ghost, noboundary},
};

fn p(e: f64) -> f64 {
    e / 3.0
}

fn dpde(_e: f64) -> f64 {
    1.0 / 3.0
}

fn constraints([t00, t01, e]: [f64; 3]) -> [f64; 5] {
    let t00 = t00.max(t01.abs());
    let e = e.max(1e-100);
    let ut = ((t00 + p(e)) / (e + p(e))).sqrt().max(1.0);
    let ux = t01 / ((e + p(e)) * ut);
    let ut = (1.0 + ux * ux).sqrt();
    [t00, t01, e, ut, ux]
}

fn eigenvalues([_t00, _t01, e, ut, ux]: [f64; 5]) -> f64 {
    let vs2 = dpde(e);
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}

pub fn f00([_, _, e, ut, _]: [f64; 5]) -> f64 {
    (e + p(e)) * ut * ut - p(e)
}
pub fn f01([_, _, e, ut, ux]: [f64; 5]) -> f64 {
    (e + p(e)) * ut * ux
}
fn f11([_, _, e, _, ux]: [f64; 5]) -> f64 {
    (e + p(e)) * ux * ux + p(e)
}

fn flux<const V: usize>(
    [_ov, v]: [&[[[f64; 3]; V]; 1]; 2],
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    er: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
    tocomp: ToCompute,
) -> [f64; 3] {
    let [t00, t01, e, _ut, _ux] = constraints(v[0][bound[0](pos[0], V)]);
    let rt0 = if tocomp.integrated() || tocomp.all() {
        let theta = 1.5;

        let divf0 = kt(
            v,
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
    let re = match tocomp {
        ToCompute::All => t00 - t01 * t01 / (t00 + p(e)),
        ToCompute::NonIntegrated => newton(er, e, |e| t00 - (t01 * t01) / (t00 + p(e)) - e),
        ToCompute::Integrated => 0.0,
    };
    [rt0[0], rt0[1], re]
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
    init: impl Fn(f64) -> [f64; 3],
) -> ([[[f64; 3]; V]; 1], f64) {
    let mut vs = [[[0.0; 3]; V]];
    let names = ["t00", "t01", "e", "ut", "ux"];
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
    let context = Context {
        fun: &flux,
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
    let (vs, t, cost, tsteps) = run(context, name, integration, &names, &constraints);
    println!("cost: {}, tsteps: {}", cost, tsteps);
    (vs, t)
}
