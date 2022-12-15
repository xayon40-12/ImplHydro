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
    let e = e.max(1e-15).min(10.0);
    let ut = ((t00 + p(e)) / (e + p(e))).sqrt();
    let ux = t01 / ((e + p(e)) * ut);
    [t00, t01, e, ut, ux]
}

fn eigenvalues([_t00, _t01, e, ut, ux]: [f64; 5]) -> f64 {
    let vs2 = dpde(e);
    let a = (ut * ux * (1.0 - vs2)).abs();
    let b = ((ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2).sqrt();
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a + b) / d
}

fn f01([_, _, e, ux, ut]: [f64; 5]) -> f64 {
    (e + p(e)) * ut * ux
}
fn f11([_, _, e, ux, _]: [f64; 5]) -> f64 {
    (e + p(e)) * ux * ux + p(e)
}

fn flux<const V: usize>(v: &[[[f64; 3]; V]; 1], bound: &[Boundary; 2], pos: [i32; 2]) -> [f64; 3] {
    let dx = 0.1;
    let theta = 1.1;

    let [t00, t01, e, _ut, _ux] = constraints(v[0][bound[0](pos[0], V)]);
    let re = t00 - t01 * t01 / (t00 + p(e));
    let rt0 = kt(
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

    [rt0[0], rt0[1], re]
}

pub fn hydro() {
    const V: usize = 10;
    let mut vs = [[[0.0; V]], [[0.0; V]], [[0.0; V]]];
    let integrated = [true, true, false];
    for i in 0..V {
        if i < V / 2 {
            vs[0][0][i] = 10.0;
            vs[2][0][i] = 10.0;
        } else {
            vs[0][0][i] = 1.0;
            vs[2][0][i] = 1.0;
        }
    }
    let context = Context {
        fun: &flux,
        boundary: &[&ghost, &noboundary], // use noboundary to emulate 1D
        local_interaction: [2, 0],        // use a distance of 0 to emulate 1D
        vs,
        integrated,
        r: [[1.0]],
        dt: 0.01,
        er: 1e-10,
        tbeg: 0.0,
        tend: 10.0,
    };
    let (vs, tot_f, tsteps) = run(context);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
    println!("t00: {:?}", vs[0][0]);
    println!("t01: {:?}", vs[1][0]);
    println!("e: {:?}", vs[2][0]);
}
