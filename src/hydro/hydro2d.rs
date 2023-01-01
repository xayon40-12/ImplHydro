use crate::solver::{
    context::{Boundary, Context, Integration, ToCompute},
    kt::{kt, Dir},
    newton::newton,
    run,
    utils::ghost,
    Constraints,
};

use super::{solve_v, Pressure};

fn gen_constraints<'a>(
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn([f64; 4]) -> [f64; 10] + 'a + Sync> {
    Box::new(|[t00, t01, t02, v]| {
        let m = (t01 * t01 + t02 * t02).sqrt();
        let t00 = t00.max(m);
        let e = (t00 - m * v).max(1e-100);
        let pe = p(e);
        let v = v.max(0.0).min(1.0);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();
        [t00, t01, t02, v, e, pe, dpde(e), ut, ux, uy]
    })
}

fn eigenvaluesx([_t00, _t01, _t02, _v, _e, _pe, dpde, ut, ux, _uy]: [f64; 10]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy([t00, t01, t02, v, e, pe, dpde, ut, ux, uy]: [f64; 10]) -> f64 {
    eigenvaluesx([t00, t01, t02, v, e, pe, dpde, ut, uy, ux])
}

pub fn f00([_, _, _, _, e, pe, _, ut, _, _]: [f64; 10]) -> f64 {
    (e + pe) * ut * ut - pe
}
pub fn f01([_, _, _, _, e, pe, _, ut, ux, _]: [f64; 10]) -> f64 {
    (e + pe) * ut * ux
}
pub fn f02([_, _, _, _, e, pe, _, ut, _, uy]: [f64; 10]) -> f64 {
    (e + pe) * uy * ut
}

fn f11([_, _, _, _, e, pe, _, _, ux, _]: [f64; 10]) -> f64 {
    (e + pe) * ux * ux + pe
}
fn f12([_, _, _, _, e, pe, _, _, ux, uy]: [f64; 10]) -> f64 {
    (e + pe) * ux * uy
}

fn f21([_, _, _, _, e, pe, _, _, ux, uy]: [f64; 10]) -> f64 {
    (e + pe) * uy * ux
}
fn f22([_, _, _, _, e, pe, _, _, _, uy]: [f64; 10]) -> f64 {
    (e + pe) * uy * uy + pe
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 4]; V]; V]; 2],
    constraints: Constraints<4, 10>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    er: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    opt: &Coordinate,
    tocomp: ToCompute,
) -> [f64; 4] {
    let [t00, t01, t02, v, e, pe, _dpde, ut, ux, uy] =
        constraints(vs[bound[1](pos[1], V)][bound[0](pos[0], V)]);
    let rt0 = if tocomp.integrated() || tocomp.all() {
        let theta = 1.5;

        let divf1 = kt(
            vs,
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
            vs,
            bound,
            pos,
            Dir::Y,
            [&f02, &f21, &f22],
            &constraints,
            &eigenvaluesy,
            dx,
            theta,
        );

        // let re = newton(1e-10, e, |e| {
        //     t00 - (t01 * t01 + t02 * t02) / (t00 + p(e)) - e
        // }); // WARNING: the energy density must be solved in the flux and not in the constraints, else there are oscillations
        let source: [f64; 3] = match opt {
            Coordinate::Cartesian => [0.0; 3],
            Coordinate::Milne => [
                (e + pe) * ut * ut / t,
                (e + pe) * ut * ux / t,
                (e + pe) * ut * uy / t,
            ],
        };
        [
            -divf1[0] - divf2[0] - source[0],
            -divf1[1] - divf2[1] - source[1],
            -divf1[2] - divf2[2] - source[2],
        ]
    } else {
        [0.0, 0.0, 0.0]
    };
    let m = (t01 * t01 + t02 * t02).sqrt();
    let sv = solve_v(t00, m);
    let rv = match tocomp {
        ToCompute::All => sv(v),
        ToCompute::NonIntegrated => newton(er, v, |v| sv(v) - v),
        ToCompute::Integrated => 0.0,
    };
    [rt0[0], rt0[1], rt0[2], rv]
}

pub fn hydro2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: ([[f64; S]; S], Option<[f64; S]>),
    opt: Coordinate,
    integration: Integration,
    p: Pressure,
    dpde: Pressure,
    init: impl Fn(f64, f64) -> [f64; 4],
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    let mut vs = [[[0.0; 4]; V]; V];
    let names = [
        "t00", "t01", "t02", "v", "e", "pe", "dpde", "ut", "ux", "uy",
    ];
    let mut k = [[[[0.0; 4]; V]; V]; S];
    let integrated = [true, true, true, false];
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init(x, y);
            for s in 0..S {
                k[s][j][i][3] = vs[j][i][3];
            }
        }
    }
    let constraints = gen_constraints(&p, &dpde);
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        local_interaction: [2, 2],   // use a distance of 0 to emulate 1D
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
        opt,
    };
    run(context, name, integration, &names, &constraints)
}
