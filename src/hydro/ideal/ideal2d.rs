use crate::{
    hydro::{
        utils::{eigenvaluesk, Coordinate},
        C_IDEAL_2D, F_IDEAL_2D, HBARC,
    },
    solver::{
        context::{Arr, BArr, Boundary, Context, DIM},
        run,
        space::{kt::kt, Eigenvalues, FluxInfo, InDir::*},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
    FLOAT,
};
use boxarray::boxarray;

use crate::hydro::{solve_v, Eos, Init2D, VOID};

pub fn init_from_entropy_density_2d<'a, const VX: usize, const VY: usize>(
    t0: FLOAT,
    s: &'a [[FLOAT; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
    entropy: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (FLOAT, FLOAT)) -> [FLOAT; F_IDEAL_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let s = s[j][i] / HBARC / t0;
        let e = newton(s, |e| entropy(e) - s, |e| e.max(0.0).min(1e10)).max(VOID);
        // let e = s;
        let vars = [e, p(e), dpde(e), 1.0, 0.0, 0.0];
        f0(t0, vars)
    })
}

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    coord: &Coordinate,
) -> Box<dyn Fn(FLOAT, [FLOAT; F_IDEAL_2D]) -> ([FLOAT; F_IDEAL_2D], [FLOAT; C_IDEAL_2D]) + 'a + Sync>
{
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
        let v = newton(0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let e = (t00 - m * v).max(VOID);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();
        ([t * t00, t * t01, t * t02], [e, pe, dpde(e), ut, ux, uy])
    })
}

fn eigenvaluesx(_t: FLOAT, [_e, _pe, dpde, ut, ux, _uy]: [FLOAT; C_IDEAL_2D]) -> FLOAT {
    eigenvaluesk(dpde, ut, ux)
}
fn eigenvaluesy(_t: FLOAT, [_e, _pe, dpde, ut, _ux, uy]: [FLOAT; C_IDEAL_2D]) -> FLOAT {
    eigenvaluesk(dpde, ut, uy)
}

pub fn f0(t: FLOAT, [e, pe, _, ut, ux, uy]: [FLOAT; C_IDEAL_2D]) -> [FLOAT; F_IDEAL_2D] {
    [
        t * ((e + pe) * ut * ut - pe),
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * uy * ut),
    ]
}

fn f1(t: FLOAT, [e, pe, _, ut, ux, uy]: [FLOAT; C_IDEAL_2D]) -> [FLOAT; F_IDEAL_2D] {
    [
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ux * ux + pe),
        t * ((e + pe) * ux * uy),
    ]
}

fn f2(t: FLOAT, [e, pe, _, ut, ux, uy]: [FLOAT; C_IDEAL_2D]) -> [FLOAT; F_IDEAL_2D] {
    [
        t * ((e + pe) * uy * ut),
        t * ((e + pe) * uy * ux),
        t * ((e + pe) * uy * uy + pe),
    ]
}

fn flux<const V: usize>(
    _k: &Arr<F_IDEAL_2D, V, V, 1>,
    [_ov, vs]: [&Arr<F_IDEAL_2D, V, V, 1>; 2],
    [_otrs, trs]: [&Arr<C_IDEAL_2D, V, V, 1>; 2],
    constraints: Constraint<F_IDEAL_2D, C_IDEAL_2D>,
    bound: Boundary<F_IDEAL_2D, V, V, 1>,
    pos: [i32; DIM],
    dxs: [FLOAT; DIM],
    [_ot, t]: [FLOAT; 2],
    [_dt, _cdt]: [FLOAT; 2],
    (coord, [eigx, eigy]): &(Coordinate, [Eigenvalues<C_IDEAL_2D>; 2]),
) -> [FLOAT; F_IDEAL_2D] {
    let theta = 1.1;

    let t = match coord {
        Coordinate::Cartesian => 1.0,
        Coordinate::Milne => t,
    };

    let pre = &|_t: FLOAT, vs: [FLOAT; F_IDEAL_2D]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        [m, t01, t02]
    };
    let post = &|_t: FLOAT, vs: [FLOAT; F_IDEAL_2D]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let t00 = (m * m + k).sqrt();
        [t00, t01, t02]
    };

    let diff = kt;
    let flux_infos = [
        X(FluxInfo {
            flux: &f1,
            secondary: &|_, _| [],
            eigenvalues: *eigx,
        }),
        Y(FluxInfo {
            flux: &f2,
            secondary: &|_, _| [],
            eigenvalues: *eigy,
        }),
    ];
    let [(dxf, _), (dyf, _)] = diff(
        (vs, trs),
        bound,
        pos,
        t,
        flux_infos,
        constraints,
        pre,
        post,
        dxs,
        theta,
    );

    let s: FLOAT = match coord {
        Coordinate::Cartesian => 0.0,
        Coordinate::Milne => {
            let y = pos[1] as usize;
            let x = pos[0] as usize;
            let [_e, pe, _dpde, _ut, _ux, _uy] = trs[0][y][x];
            pe
        }
    };
    [-dxf[0] - dyf[0] - s, -dxf[1] - dyf[1], -dxf[2] - dyf[2]]
}

pub fn momentum_anisotropy<const VX: usize, const VY: usize>(
    t: FLOAT,
    _vs: &Arr<F_IDEAL_2D, VX, VY, 1>,
    tran: &Arr<C_IDEAL_2D, VX, VY, 1>,
) -> Vec<FLOAT> {
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
    name: &(&str, usize),
    maxdt: FLOAT,
    t: FLOAT,
    tend: FLOAT,
    dx: FLOAT,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    init: Init2D<F_IDEAL_2D>,
    save_raw: Option<FLOAT>,
) -> Option<(
    (BArr<F_IDEAL_2D, V, V, 1>, BArr<C_IDEAL_2D, V, V, 1>),
    FLOAT,
    usize,
    usize,
)> {
    let coord = Coordinate::Milne;
    let constraints = gen_constraints(&p, &dpde, &coord);
    let mut vs: Box<[[[[FLOAT; F_IDEAL_2D]; V]; V]; 1]> = boxarray(0.0);
    let mut trs: Box<[[[[FLOAT; C_IDEAL_2D]; V]; V]; 1]> = boxarray(0.0);

    let names = (["t00", "t01", "t02"], ["e", "pe", "dpde", "ut", "ux", "uy"]);
    let k: Box<[[[[[FLOAT; F_IDEAL_2D]; V]; V]; 1]; S]> = boxarray(0.0);
    let v2 = ((V - 1) as FLOAT) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as FLOAT - v2) * dx;
            let y = (j as FLOAT - v2) * dx;
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
        dxs: [dx, dx, 0.0],
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
    let err_thr = |_t: FLOAT,
                   vs: &[[[[FLOAT; F_IDEAL_2D]; V]; V]; 1],
                   _trs: &[[[[FLOAT; C_IDEAL_2D]; V]; V]; 1]| {
        let m = vs[0]
            .iter()
            .flat_map(|v| v.iter().map(|v| v[0]))
            .sum::<FLOAT>()
            / (V * V) as FLOAT;
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
        save_raw,
    )
}
