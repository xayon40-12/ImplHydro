use crate::solver::{
    context::{Boundary, Context, Integration},
    run,
    space::{kt::kt, Dir, Eigenvalues},
    time::{newton::newton, schemes::Scheme},
    utils::{ghost, zeros},
    Constraint, Observable,
};

use crate::hydro::{solve_v, Eos, Init2D, VOID};

fn gen_constraints<'a>(
    er: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
    coord: &Coordinate,
) -> Box<dyn Fn(f64, [f64; 3]) -> ([f64; 3], [f64; 6]) + 'a + Sync> {
    let st = match coord {
        Coordinate::Cartesian => |_| 1.0,
        Coordinate::Milne => |t| t,
    };
    Box::new(move |t, [t00, t01, t02]| {
        let t = st(t);
        let t00 = t00 / t;
        let t01 = t01 / t;
        let t02 = t02 / t;
        let m = (t01 * t01 + t02 * t02).sqrt();
        let t00 = t00.max(m * (1.0 + 1e-15));
        let sv = solve_v(t00, m, p);
        let v = newton(er, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let e = (t00 - m * v).max(VOID);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();
        ([t * t00, t * t01, t * t02], [e, pe, dpde(e), ut, ux, uy])
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

pub fn f0(t: f64, [e, pe, _, ut, ux, uy]: [f64; 6]) -> [f64; 3] {
    [
        t * ((e + pe) * ut * ut - pe),
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * uy * ut),
    ]
}

fn f1(t: f64, [e, pe, _, ut, ux, uy]: [f64; 6]) -> [f64; 3] {
    [
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ux * ux + pe),
        t * ((e + pe) * ux * uy),
    ]
}

fn f2(t: f64, [e, pe, _, ut, ux, uy]: [f64; 6]) -> [f64; 3] {
    [
        t * ((e + pe) * uy * ut),
        t * ((e + pe) * uy * ux),
        t * ((e + pe) * uy * uy + pe),
    ]
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 3]; V]; V]; 2],
    [_otrs, trs]: [&[[[f64; 6]; V]; V]; 2],
    constraints: Constraint<3, 6>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    (coord, [eigx, eigy]): &(Coordinate, [Eigenvalues<6>; 2]),
) -> [f64; 3] {
    let theta = 1.1;
    // let theta = 2.0;

    let t = match coord {
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

    // let pre = &id_flux_limiter;
    // let post = &id_flux_limiter;

    let diff = kt;
    let (divf1, _) = diff(
        vs,
        bound,
        pos,
        Dir::X,
        t,
        &f1,
        (&|_, _| [], []),
        constraints,
        *eigx,
        pre,
        post,
        dx,
        theta,
    );
    let (divf2, _) = diff(
        vs,
        bound,
        pos,
        Dir::Y,
        t,
        &f2,
        (&|_, _| [], []),
        constraints,
        *eigy,
        pre,
        post,
        dx,
        theta,
    );

    let s: f64 = match coord {
        Coordinate::Cartesian => 0.0,
        Coordinate::Milne => {
            let y = bound[1](pos[1], V);
            let x = bound[0](pos[0], V);
            let [_e, pe, _dpde, _ut, _ux, _uy] = trs[y][x];
            pe
        }
    };
    [
        -divf1[0] - divf2[0] - s,
        -divf1[1] - divf2[1],
        -divf1[2] - divf2[2],
    ]
}

pub fn momentum_anysotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &[[[f64; 3]; VX]; VY],
    tran: &[[[f64; 6]; VX]; VY],
) -> Vec<f64> {
    let mut mt11 = 0.0;
    let mut mt22 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            mt11 += f1(t, tran[j][i])[1];
            mt22 += f2(t, tran[j][i])[2];
        }
    }
    let anysotropy = (mt11 - mt22) / (mt11 + mt22);
    vec![anysotropy]
}

pub fn ideal2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    init: Init2D<3>,
) -> Option<(([[[f64; 3]; V]; V], [[[f64; 6]; V]; V]), f64, usize, usize)> {
    let coord = Coordinate::Milne;
    let constraints = gen_constraints(er, &p, &dpde, &coord);
    let schemename = r.name;
    let mut vs = [[[0.0; 3]; V]; V];
    let mut trs = [[[0.0; 6]; V]; V];

    let names = (["t00", "t01", "t02"], ["e", "pe", "dpde", "ut", "ux", "uy"]);
    let mut k = [[[[0.0; 3]; V]; V]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            (vs[j][i], trs[j][i]) = constraints(t, init((i, j), (x, y)));
        }
    }
    match r.integration {
        Integration::Explicit => k[S - 1] = vs, // prepare k[S-1] so that it can be use as older time for time derivatives (in this case approximate time derivatives to be zero at initial time)
        Integration::FixPoint(_) => {}
    }
    let integration = r.integration;

    let eigx = Eigenvalues::Analytical(&eigenvaluesx);
    let eigy = Eigenvalues::Analytical(&eigenvaluesy);

    // let eigconstr = |_t: f64, [_, _, dpde, _, _, _]: [f64; 6], eig: f64| eig.max(dpde).min(1.0);
    // let eigx: Eigenvalues<6> = Eigenvalues::ApproxConstraint(&eigconstr);
    // let eigy = eigx;
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        post_constraints: None,
        local_interaction: [1, 1], // use a distance of 0 to emulate 1D
        vstrs: (vs, trs),
        total_diff_vs: zeros(&vs),
        k,
        r,
        dt: 1e10,
        dx,
        maxdt,
        er,
        t,
        ot: t,
        t0: t,
        tend,
        opt: (coord, [eigx, eigy]),
        p,
        dpde,
    };

    let observables: [Observable<3, 6, V, V>; 1] =
        [("momentum_anysotropy", &momentum_anysotropy::<V, V>)];

    run(
        context,
        name,
        &schemename,
        integration,
        &names,
        &observables,
    )
}
