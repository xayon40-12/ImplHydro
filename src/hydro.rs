use crate::{
    newton::{Boundary, Context},
    solver::run,
    utils::{flux_limiter, ghost, noboundary, periodic},
};

fn p(e: f64) -> f64 {
    e / 3.0
}

fn constraints([t00, t01, e]: [f64; 3]) -> [f64; 5] {
    let e = e.max(1e-15).min(10.0);
    let ut = ((t00 + p(e)) / (e + p(e))).sqrt();
    let ux = t01 / ((e + p(e)) * ut);
    [t00, t01, e, ut, ux]
}

fn f01([_, _, e, ux, ut]: [f64; 5]) -> f64 {
    (e + p(e)) * ut * ux
}
fn f11([_, _, e, ux, _]: [f64; 5]) -> f64 {
    (e + p(e)) * ux * ux + p(e)
}

fn flux<const V: usize>(
    v: &[[[f64; 3]; V]; 1],
    [bound, _]: &[Boundary; 2],
    [i, _]: [i32; 2],
) -> [f64; 3] {
    let ivdx = 10.0f64;
    let theta = 1.1;

    let vals: Vec<[f64; 5]> = (-2..=2)
        .map(|l| constraints(v[0][bound(i - l, V)]))
        .collect();
    let j01: Vec<f64> = vals.iter().map(|v| f01(*v)).collect();
    let j11: Vec<f64> = vals.iter().map(|v| f11(*v)).collect();

    let rt00 = ivdx * flux_limiter(theta, j01[1], j01[2], j01[3]);
    let rt01 = ivdx * flux_limiter(theta, j11[1], j11[2], j11[3]);

    let t00 = vals[2][0];
    let t01 = vals[2][1];
    let e = vals[2][2];
    let re = t00 - t01 * t01 / (t00 + p(e));

    [rt00, rt01, re]
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
