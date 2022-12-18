use crate::{
    kt::{kt, Dir},
    newton::{Boundary, Context},
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
    let e = e.max(1e-15).min(100.0);
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

fn f01([_, _, e, ut, ux]: [f64; 5]) -> f64 {
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
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
) -> [f64; 3] {
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

    let [t00, t01, e, _ut, _ux] = constraints(v[0][bound[0](pos[0], V)]);
    let re = t00 - t01 * t01 / (t00 + p(e));
    let rt0 = [-divf0[0], -divf0[1]];
    [rt0[0], rt0[1], re]
}

pub fn hydro1d() {
    const V: usize = 50;
    let mut vs = [[[0.0; 3]; V]];
    let names = ["t00", "t01", "e", "ut", "ux"];
    let k = [[[[0.0; 3]; V]]];
    let integrated = [true, true, false];
    for i in 0..V {
        let e = if i == V / 2 { 10.0 } else { 1.0 };
        vs[0][i] = [e, 0.0, e];
    }
    let context = Context {
        fun: &flux,
        boundary: &[&ghost, &noboundary], // use noboundary to emulate 1D
        local_interaction: [2, 0],        // use a distance of 0 to emulate 1D
        vs,
        k,
        integrated,
        r: [[1.0]],
        dt: 0.1,
        dx: 0.2,
        maxdt: 0.1,
        er: 1e-5,
        t: 0.0,
        tend: 5.0,
        opt: (),
    };
    let (_vs, cost, tsteps) = run(context, &names, &constraints);
    println!("cost: {}, tsteps: {}", cost, tsteps);
}
