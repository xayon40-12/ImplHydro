use crate::{
    boxarray::boxarray,
    hydro::{utils::eigenvaluesk, C_IDEAL_1D, F_IDEAL_1D},
    solver::{
        context::{Arr, BArr, Boundary, Context, DIM},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint,
    },
};

use crate::hydro::{solve_v, Eos, Init1D, VOID};

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn(f64, [f64; F_IDEAL_1D]) -> ([f64; F_IDEAL_1D], [f64; C_IDEAL_1D]) + 'a + Sync> {
    Box::new(move |_t, [t00, t01]| {
        let m = t01.abs();
        let t00 = t00.max(m * (1.0 + 1e-15));
        let sv = solve_v(t00, m, p);
        let v = newton(1e-10, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let e = (t00 - m * v).max(VOID);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux).sqrt();
        ([t00, t01], [e, pe, dpde(e), ut, ux])
    })
}

fn eigenvalues(_t: f64, [_e, _pe, dpde, ut, ux]: [f64; C_IDEAL_1D]) -> f64 {
    eigenvaluesk(dpde, ut, ux)
}

pub fn f0(_t: f64, [e, pe, _, ut, ux]: [f64; C_IDEAL_1D]) -> [f64; F_IDEAL_1D] {
    [(e + pe) * ut * ut - pe, (e + pe) * ut * ux]
}

fn f1(_t: f64, [e, pe, _, ut, ux]: [f64; C_IDEAL_1D]) -> [f64; F_IDEAL_1D] {
    [(e + pe) * ut * ux, (e + pe) * ux * ux + pe]
}

fn flux<const V: usize>(
    [_ov, vs]: [&Arr<F_IDEAL_1D, V, 1, 1>; 2],
    [_otrs, trs]: [&Arr<C_IDEAL_1D, V, 1, 1>; 2],
    constraints: Constraint<F_IDEAL_1D, C_IDEAL_1D>,
    bound: Boundary<F_IDEAL_1D, V, 1, 1>,
    pos: [i32; DIM],
    dx: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
) -> [f64; 2] {
    let theta = 1.1;

    let pre = &|_t: f64, vs: [f64; 2]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let k = t01 * t01;
        let m = (t00 * t00 - k).sqrt();
        [m, t01]
    };
    let post = &|_t: f64, vs: [f64; 2]| {
        let m = vs[0];
        let t01 = vs[1];
        let k = t01 * t01;
        let t00 = (m * m + k).sqrt();
        [t00, t01]
    };

    let diff = kt;

    let (divf0, _) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::X,
        1.0,
        &f1,
        &|_, _| [],
        constraints,
        Eigenvalues::Analytical(&eigenvalues),
        pre,
        post,
        dx,
        theta,
    );

    [-divf0[0], -divf0[1]]
}

pub fn ideal1d<const V: usize, const S: usize>(
    name: &(&str, usize),
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    init: Init1D<2>,
    save_raw: bool,
) -> Option<(
    (BArr<F_IDEAL_1D, V, 1, 1>, BArr<C_IDEAL_1D, V, 1, 1>),
    f64,
    usize,
    usize,
)> {
    let constraints = gen_constraints(p, dpde);
    let mut vs: Box<[[[[f64; F_IDEAL_1D]; V]; 1]; 1]> = boxarray(0.0);
    let mut trs: Box<[[[[f64; C_IDEAL_1D]; V]; 1]; 1]> = boxarray(0.0);
    let names = (["t00", "t01"], ["e", "pe", "dpde", "ut", "ux"]);
    let k: Box<[[[[[f64; 2]; V]; 1]; 1]; S]> = boxarray(0.0);
    let v2 = ((V - 1) as f64) / 2.0;
    for i in 0..V {
        let x = (i as f64 - v2) * dx;
        vs[0][0][i] = init(i, x);
        (vs[0][0][i], trs[0][0][i]) = constraints(t, vs[0][0][i]);
    }
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &ghost, // use noboundary to emulate 1D
        post_constraints: None,
        local_interaction: [1, 0, 0], // use a distance of 0 to emulate 1D
        ovstrs: (vs.clone(), trs.clone()),
        vstrs: (vs, trs),
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
        opt: (),
        p,
        dpde,
        freezeout_energy: None,
    };

    let e = 1e-3;
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_IDEAL_1D]; V]; 1]; 1],
                   _trs: &[[[[f64; C_IDEAL_1D]; V]; 1]; 1]| {
        let m = vs[0][0].iter().map(|v| v[0]).sum::<f64>() / V as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order)
    };

    run(
        context,
        name,
        crate::hydro::Viscosity::Ideal,
        &names,
        &[],
        &err_thr,
        save_raw,
    )
}
