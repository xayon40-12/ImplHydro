use crate::{
    hydro::{Viscosity, C_BOTH_2D, F_BOTH_2D, HBARC},
    solver::{
        context::{Boundary, Context, Integration},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{solve_v, Eos, Init2D, VOID};

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    _temp: Eos<'a>,
    _implicit: bool,
) -> Box<dyn Fn(f64, [f64; F_BOTH_2D]) -> ([f64; F_BOTH_2D], [f64; C_BOTH_2D]) + 'a + Sync> {
    Box::new(move |t, [t00, t01, t02, _, _, _, _]| {
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
        (
            [t * t00, t * t01, t * t02, 0.0, 0.0, 0.0, 0.0],
            [
                e,
                pe,
                dpde(e),
                ut,
                ux,
                uy,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )
    })
}

fn eigenvaluesx(_t: f64, [_, _, dpde, ut, ux, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_2D]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy(
    t: f64,
    [e, pe, dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33, ppi]: [f64; C_BOTH_2D],
) -> f64 {
    eigenvaluesx(
        t,
        [
            e, pe, dpde, ut, uy, ux, pi00, pi01, pi02, pi11, pi12, pi22, pi33, ppi,
        ],
    )
}

pub fn fitutpi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    [
        t * ((e + pe) * ut * ut - pe), // no pi in Ttt because it is substracted in the flux
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * uy * ut),
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}
fn f1(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    [
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ux * ux + pe),
        t * ((e + pe) * ux * uy),
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}
fn f2(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    [
        t * ((e + pe) * uy * ut),
        t * ((e + pe) * uy * ux),
        t * ((e + pe) * uy * uy + pe),
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; F_BOTH_2D]; V]; V]; 2],
    [_, trs]: [&[[[f64; C_BOTH_2D]; V]; V]; 2],
    constraints: Constraint<F_BOTH_2D, C_BOTH_2D>,
    bound: Boundary<F_BOTH_2D, V, V>,
    pos: [i32; 2],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _: &((f64, f64, f64), (f64, f64, f64), Eos, f64),
) -> [f64; F_BOTH_2D] {
    let theta = 1.1;

    let pre = &|_t: f64, vs: [f64; F_BOTH_2D]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        [m, t01, t02, 0.0, 0.0, 0.0, 0.0]
    };
    let post = &|_t: f64, vs: [f64; F_BOTH_2D]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let t00 = (m * m + k).sqrt();
        [t00, t01, t02, 0.0, 0.0, 0.0, 0.0]
    };

    // let pre = &id_flux_limiter;
    // let post = &id_flux_limiter;

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
        Eigenvalues::Analytical(&eigenvaluesx),
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
        Eigenvalues::Analytical(&eigenvaluesy),
        pre,
        post,
        dx,
        theta,
    );

    let s: f64 = {
        let y = pos[1] as usize;
        let x = pos[0] as usize;
        let [_e, pe, _dpde, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi] =
            trs[y][x];
        pe
    };
    [
        -divf1[0] - divf2[0] - s,
        -divf1[1] - divf2[1],
        -divf1[2] - divf2[2],
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

pub fn momentum_anysotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &[[[f64; F_BOTH_2D]; VX]; VY],
    tran: &[[[f64; C_BOTH_2D]; VX]; VY],
) -> Vec<f64> {
    let mut mt11 = 0.0;
    let mut mt22 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            mt11 += f1(t, tran[j][i])[1];
            mt22 += f2(t, tran[j][i])[2];
        }
    }
    let anysotropy = if mt11 + mt22 == 0.0 {
        0.0
    } else {
        (mt11 - mt22) / (mt11 + mt22)
    };
    vec![anysotropy]
}

// viscous hydro is in Milne coordinates
pub fn viscous2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    temperature: Eos,
    init: Init2D<F_BOTH_2D>,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    shear_temp_cut: f64,
    freezeout_temp: f64,
) -> Option<(
    ([[[f64; F_BOTH_2D]; V]; V], [[[f64; C_BOTH_2D]; V]; V]),
    f64,
    usize,
    usize,
)> {
    let implicit = match r.integration {
        Integration::FixPoint => true,
        Integration::Explicit => false,
    };
    let constraints = gen_constraints(&p, &dpde, temperature, implicit);

    let mut vs = [[[0.0; F_BOTH_2D]; V]; V];
    let mut trs = [[[0.0; C_BOTH_2D]; V]; V];
    let names = (
        ["tt00", "tt01", "tt02", "utpi11", "utpi12", "utpi22", "utPi"],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "pi00", "pi01", "pi02", "pi11", "pi12", "pi22",
            "pi33", "Pi",
        ],
    );
    let k = [[[[0.0; F_BOTH_2D]; V]; V]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    let mut max_e = 0.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init((i, j), (x, y));
            trs[j][i][0] = vs[j][i][0];
            (vs[j][i], trs[j][i]) = constraints(t, vs[j][i]);
            max_e = trs[j][i][0].max(max_e);
        }
    }

    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &ghost, // TODO use better boundary
        post_constraints: None,
        local_interaction: [1, 1], // use a distance of 0 to emulate 1D
        vstrs: (vs, trs),
        ovstrs: (vs, trs),
        total_diff_vs: zeros(&vs),
        k,
        r,
        dt: 1e10,
        dx,
        maxdt,
        t,
        ot: t - 1.0,
        t0: t,
        tend,
        opt: (etaovers, zetaovers, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: None,
    };

    let e = 2e-4;
    let err_thr = |_t: f64, _vs: &[[[f64; F_BOTH_2D]; V]; V], _trs: &[[[f64; C_BOTH_2D]; V]; V]| {
        let m = _vs.iter().flat_map(|v| v.iter().map(|v| v[0])).sum::<f64>() / (V * V) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order)
    };

    let observables: [Observable<F_BOTH_2D, C_BOTH_2D, V, V>; 1] =
        [("momentum_anysotropy", &momentum_anysotropy::<V, V>)];

    run(
        context,
        name,
        crate::hydro::Viscosity::Ideal,
        &names,
        &observables,
        &err_thr,
    )
}
