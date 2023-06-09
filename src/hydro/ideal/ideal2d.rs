use crate::{
    boxarray,
    hydro::{C_IDEAL_2D, F_IDEAL_2D},
    solver::{
        context::{Arr, BArr, Boundary, Context, DIM},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{solve_v, Eos, Init2D, VOID};

pub fn init_from_entropy_density_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    s: &'a [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_IDEAL_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let s = s[j][i].max(VOID);
        let e = s;
        let vars = [e, p(e), dpde(e), 1.0, 0.0, 0.0];
        f0(t0, vars)
    })
}

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    coord: &Coordinate,
) -> Box<dyn Fn(f64, [f64; F_IDEAL_2D]) -> ([f64; F_IDEAL_2D], [f64; C_IDEAL_2D]) + 'a + Sync> {
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
        let v = newton(1e-10, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let e = (t00 - m * v).max(VOID);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();
        ([t * t00, t * t01, t * t02], [e, pe, dpde(e), ut, ux, uy])
    })
}

fn eigenvaluesx(_t: f64, [_e, _pe, dpde, ut, ux, _uy]: [f64; C_IDEAL_2D]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy(t: f64, [e, pe, dpde, ut, ux, uy]: [f64; C_IDEAL_2D]) -> f64 {
    eigenvaluesx(t, [e, pe, dpde, ut, uy, ux])
}

pub fn f0(t: f64, [e, pe, _, ut, ux, uy]: [f64; C_IDEAL_2D]) -> [f64; F_IDEAL_2D] {
    [
        t * ((e + pe) * ut * ut - pe),
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * uy * ut),
    ]
}

fn f1(t: f64, [e, pe, _, ut, ux, uy]: [f64; C_IDEAL_2D]) -> [f64; F_IDEAL_2D] {
    [
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ux * ux + pe),
        t * ((e + pe) * ux * uy),
    ]
}

fn f2(t: f64, [e, pe, _, ut, ux, uy]: [f64; C_IDEAL_2D]) -> [f64; F_IDEAL_2D] {
    [
        t * ((e + pe) * uy * ut),
        t * ((e + pe) * uy * ux),
        t * ((e + pe) * uy * uy + pe),
    ]
}

#[derive(Clone)]
pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&Arr<F_IDEAL_2D, V, V, 1>; 2],
    [_otrs, trs]: [&Arr<C_IDEAL_2D, V, V, 1>; 2],
    constraints: Constraint<F_IDEAL_2D, C_IDEAL_2D>,
    bound: Boundary<F_IDEAL_2D, V, V, 1>,
    pos: [i32; DIM],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    (coord, [eigx, eigy]): &(Coordinate, [Eigenvalues<C_IDEAL_2D>; 2]),
) -> [f64; F_IDEAL_2D] {
    let theta = 1.1;

    let t = match coord {
        Coordinate::Cartesian => 1.0,
        Coordinate::Milne => t,
    };

    let pre = &|_t: f64, vs: [f64; F_IDEAL_2D]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        [m, t01, t02]
    };
    let post = &|_t: f64, vs: [f64; F_IDEAL_2D]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let t00 = (m * m + k).sqrt();
        [t00, t01, t02]
    };

    let diff = kt;
    let (divf1, _) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::X,
        t,
        &f1,
        &|_, _| [],
        constraints,
        *eigx,
        pre,
        post,
        dx,
        theta,
    );
    let (divf2, _) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::Y,
        t,
        &f2,
        &|_, _| [],
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
            let y = pos[1] as usize;
            let x = pos[0] as usize;
            let [_e, pe, _dpde, _ut, _ux, _uy] = trs[0][y][x];
            pe
        }
    };
    [
        -divf1[0] - divf2[0] - s,
        -divf1[1] - divf2[1],
        -divf1[2] - divf2[2],
    ]
}

pub fn momentum_anisotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &Arr<F_IDEAL_2D, VX, VY, 1>,
    tran: &Arr<C_IDEAL_2D, VX, VY, 1>,
) -> Vec<f64> {
    let mut mt11 = 0.0;
    let mut mt12 = 0.0;
    let mut mt22 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            let [_, t11, t12] = f1(t, tran[0][j][i]);
            let [_, _, t22] = f2(t, tran[0][j][i]);
            mt11 += t11;
            mt12 += t12;
            mt22 += t22;
        }
    }
    let anisotropy = if mt11 + mt22 == 0.0 {
        0.0
    } else {
        (mt11 - mt22).hypot(2.0 * mt12) / (mt11 + mt22)
    };
    vec![anisotropy]
}

pub fn ideal2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    init: Init2D<F_IDEAL_2D>,
) -> Option<(
    (BArr<F_IDEAL_2D, V, V, 1>, BArr<C_IDEAL_2D, V, V, 1>),
    f64,
    usize,
    usize,
)> {
    let coord = Coordinate::Milne;
    let constraints = gen_constraints(&p, &dpde, &coord);
    let mut vs: Box<[[[[f64; F_IDEAL_2D]; V]; V]; 1]> = boxarray(0.0);
    let mut trs: Box<[[[[f64; C_IDEAL_2D]; V]; V]; 1]> = boxarray(0.0);

    let names = (["t00", "t01", "t02"], ["e", "pe", "dpde", "ut", "ux", "uy"]);
    let k: Box<[[[[[f64; F_IDEAL_2D]; V]; V]; 1]; S]> = boxarray(0.0);
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[0][j][i] = init((i, j), (x, y));
            (vs[0][j][i], trs[0][j][i]) = constraints(t, vs[0][j][i]);
        }
    }

    let eigx = Eigenvalues::Analytical(&eigenvaluesx);
    let eigy = Eigenvalues::Analytical(&eigenvaluesy);

    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &ghost, // use noboundary to emulate 1D
        post_constraints: None,
        local_interaction: [1, 1, 0], // use a distance of 0 to emulate 1D
        vstrs: (vs.clone(), trs.clone()),
        ovstrs: (vs, trs),
        total_diff_vs: zeros(),
        k,
        r,
        dt: 1e10,
        dx,
        maxdt,
        t,
        ot: t - 1.0,
        t0: t,
        tend,
        opt: (coord, [eigx, eigy]),
        p,
        dpde,
        freezeout_energy: None,
    };

    let e = 2e-4;
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_IDEAL_2D]; V]; V]; 1],
                   _trs: &[[[[f64; C_IDEAL_2D]; V]; V]; 1]| {
        let m = vs[0]
            .iter()
            .flat_map(|v| v.iter().map(|v| v[0]))
            .sum::<f64>()
            / (V * V) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order)
    };

    let observables: [Observable<F_IDEAL_2D, C_IDEAL_2D, V, V, 1>; 1] =
        [("momentum_anisotropy", &momentum_anisotropy::<V, V>)];

    run(
        context,
        name,
        crate::hydro::Viscosity::Ideal,
        &names,
        &observables,
        &err_thr,
    )
}
