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
    let e = e.max(1e-15).min(100.0);
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

fn f10([_, _, _, e, ut, ux, _]: [f64; 7]) -> f64 {
    (e + p(e)) * ut * ux
}
fn f11([_, _, _, e, _, ux, _]: [f64; 7]) -> f64 {
    (e + p(e)) * ux * ux + p(e)
}
fn f12([_, _, _, e, _, ux, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * ux * uy
}

fn f20([_, _, _, e, ut, _, uy]: [f64; 7]) -> f64 {
    (e + p(e)) * uy * ut
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
        [&f10, &f11, &f12],
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
        [&f20, &f21, &f22],
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

pub fn hydro2d(maxdt: f64, er: f64, t: f64, tend: f64, opt: Coordinate) {
    const V: usize = 100;
    let mut vs = [[[0.0; 4]; V]; V];
    let k = [[[[0.0; 4]; V]; V]];
    let integrated = [true, true, true, false];
    for i in 0..V {
        for j in 0..V {
            let e = if i == V / 2 && j == V / 2 {
                10.0
            } else {
                1e-15
            };
            vs[i][j] = [e, 0.0, 0.0, e];
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
        dx: 0.1,
        maxdt,
        er,
        t,
        tend,
        opt,
    };
    let (vals, _tot_f, _tsteps) = run(context);
    // let tot_f = tot_f * (2 * 2 + 2 * 2 + 1 + 1);
    // println!(
    //     "tot_f: {}, tsteps: {}, ratio: {}",
    //     tot_f,
    //     tsteps,
    //     tot_f as f64 / tsteps as f64
    // );
    // let mut es = [(0.0, 0.0); V];
    // let mut vs = [(0.0, 0.0); V];
    for j in 0..V {
        for i in 0..V {
            let [_t00, _t01, _t02, e, _ut, _ux, _uy] = constraints(vals[j][i]);
            let y = (j as f32 / V as f32) * 2.0 - 1.0;
            let x = (i as f32 / V as f32) * 2.0 - 1.0;
            println!("{} {} {}", x, y, e);
            // es[v] = (x, e as _);
            // vs[v] = (x, ((ux * ux + uy * uy).sqrt() / ut) as _);
        }
        println!("");
    }
    // Chart::new(400, 200, -1.0, 1.0)
    //     .lineplot(&Shape::Lines(&es))
    //     .display();
}
