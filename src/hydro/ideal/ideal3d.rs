use crate::{
    boxarray,
    hydro::{
        utils::{eigenvaluesk, Coordinate},
        C_IDEAL_3D, F_IDEAL_3D,
    },
    solver::{
        context::{Arr, BArr, Boundary, Context, DIM},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{solve_v, Eos, Init3D, VOID};

pub fn init_from_entropy_density_3d<'a, const VX: usize, const VY: usize, const VZ: usize>(
    t0: f64,
    s: &'a [[[f64; VX]; VY]; VZ],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize, usize), (f64, f64, f64)) -> [f64; F_IDEAL_3D] + 'a> {
    Box::new(move |(i, j, k), _| {
        let s = s[k][j][i].max(VOID);
        let e = s;
        let vars = [e, p(e), dpde(e), 1.0, 0.0, 0.0, 0.0];
        f0(t0, vars)
    })
}

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    coord: &Coordinate,
) -> Box<dyn Fn(f64, [f64; F_IDEAL_3D]) -> ([f64; F_IDEAL_3D], [f64; C_IDEAL_3D]) + 'a + Sync> {
    let st = match coord {
        Coordinate::Cartesian => |_| 1.0,
        Coordinate::Milne => |t| t,
    };
    Box::new(move |t, [t00, t01, t02, t03]| {
        let t = st(t);
        let t00 = t00 / t;
        let t01 = t01 / t;
        let t02 = t02 / t;
        let t03 = t03 / t;
        let m = (t01 * t01 + t02 * t02 + t03 * t03).sqrt();
        let t00 = t00.max(m * (1.0 + 1e-15));
        let sv = solve_v(t00, m, p);
        let v = newton(1e-10, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let e = (t00 - m * v).max(VOID);
        let pe = p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let uz = t03 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy + uz * uz).sqrt();
        (
            [t * t00, t * t01, t * t02, t * t03],
            [e, pe, dpde(e), ut, ux, uy, uz],
        )
    })
}

fn eigenvaluesx(_t: f64, [_e, _pe, dpde, ut, ux, _uy, _uz]: [f64; C_IDEAL_3D]) -> f64 {
    eigenvaluesk(dpde, ut, ux)
}
fn eigenvaluesy(_t: f64, [_e, _pe, dpde, ut, _ux, uy, _uz]: [f64; C_IDEAL_3D]) -> f64 {
    eigenvaluesk(dpde, ut, uy)
}
fn eigenvaluesz(t: f64, [_e, _pe, dpde, ut, _ux, _uy, uz]: [f64; C_IDEAL_3D]) -> f64 {
    eigenvaluesk(dpde, ut, uz) / t
}

pub fn f0(t: f64, [e, pe, _, ut, ux, uy, uz]: [f64; C_IDEAL_3D]) -> [f64; F_IDEAL_3D] {
    [
        t * ((e + pe) * ut * ut - pe),
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ut * uy),
        t * ((e + pe) * ut * uz),
    ]
}

fn f1(t: f64, [e, pe, _, ut, ux, uy, uz]: [f64; C_IDEAL_3D]) -> [f64; F_IDEAL_3D] {
    [
        t * ((e + pe) * ux * ut),
        t * ((e + pe) * ux * ux + pe),
        t * ((e + pe) * ux * uy),
        t * ((e + pe) * ux * uz),
    ]
}

fn f2(t: f64, [e, pe, _, ut, ux, uy, uz]: [f64; C_IDEAL_3D]) -> [f64; F_IDEAL_3D] {
    [
        t * ((e + pe) * uy * ut),
        t * ((e + pe) * uy * ux),
        t * ((e + pe) * uy * uy + pe),
        t * ((e + pe) * uy * uz),
    ]
}

fn f3(_t: f64, [e, pe, _, ut, ux, uy, uz]: [f64; C_IDEAL_3D]) -> [f64; F_IDEAL_3D] {
    [
        (e + pe) * uz * ut,
        (e + pe) * uz * ux,
        (e + pe) * uz * uy,
        (e + pe) * uz * uz + pe,
    ]
}

fn flux<const XY: usize, const Z: usize>(
    [_ov, vs]: [&Arr<F_IDEAL_3D, XY, XY, Z>; 2],
    [_otrs, trs]: [&Arr<C_IDEAL_3D, XY, XY, Z>; 2],
    constraints: Constraint<F_IDEAL_3D, C_IDEAL_3D>,
    bound: Boundary<F_IDEAL_3D, XY, XY, Z>,
    pos: [i32; DIM],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    (coord, [eigx, eigy, eigz]): &(Coordinate, [Eigenvalues<C_IDEAL_3D>; 3]),
) -> [f64; F_IDEAL_3D] {
    let theta = 1.1;

    let t = match coord {
        Coordinate::Cartesian => 1.0,
        Coordinate::Milne => t,
    };

    let pre = &|_t: f64, vs: [f64; F_IDEAL_3D]| {
        let [t00, t01, t02, t03] = vs;
        let k = t01 * t01 + t02 * t02 + t03 * t03;
        let m = (t00 * t00 - k).sqrt();
        [m, t01, t02, t03]
    };
    let post = &|_t: f64, vs: [f64; F_IDEAL_3D]| {
        let [m, t01, t02, t03] = vs;
        let k = t01 * t01 + t02 * t02 + t03 * t03;
        let t00 = (m * m + k).sqrt();
        [t00, t01, t02, t03]
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
    let (divf3, _) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::Z,
        t,
        &f3,
        &|_, _| [],
        constraints,
        *eigz,
        pre,
        post,
        dx,
        theta,
    );

    let (tee, tte): (f64, f64) = match coord {
        Coordinate::Cartesian => (0.0, 0.0),
        Coordinate::Milne => {
            let y = pos[1] as usize;
            let x = pos[0] as usize;
            let [e, pe, _dpde, ut, _ux, _uy, uz] = trs[0][y][x];
            ((e + pe) * uz * uz + pe, (e + pe) * ut * uz)
        }
    };
    [
        -divf1[0] - divf2[0] - divf3[0] - tee,
        -divf1[1] - divf2[1] - divf3[1],
        -divf1[2] - divf2[2] - divf3[2],
        -divf1[3] - divf2[3] - divf3[3] - tte,
    ]
}

pub fn ideal3d<const XY: usize, const Z: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    init: Init3D<F_IDEAL_3D>,
) -> Option<(
    (BArr<F_IDEAL_3D, XY, XY, Z>, BArr<C_IDEAL_3D, XY, XY, Z>),
    f64,
    usize,
    usize,
)> {
    let coord = Coordinate::Milne;
    let constraints = gen_constraints(&p, &dpde, &coord);
    let mut vs: Box<[[[[f64; F_IDEAL_3D]; XY]; XY]; Z]> = boxarray(0.0f64);
    let mut trs: Box<[[[[f64; C_IDEAL_3D]; XY]; XY]; Z]> = boxarray(0.0f64);

    let names = (
        ["t00", "t01", "t02", "t03"],
        ["e", "pe", "dpde", "ut", "ux", "uy", "uz"],
    );
    let k: Box<[[[[[f64; F_IDEAL_3D]; XY]; XY]; Z]; S]> = boxarray(0.0f64);
    let v2 = ((XY - 1) as f64) / 2.0;
    let v2z = ((Z - 1) as f64) / 2.0;
    for k in 0..Z {
        for j in 0..XY {
            for i in 0..XY {
                let x = (i as f64 - v2) * dx;
                let y = (j as f64 - v2) * dx;
                let z = (k as f64 - v2z) * dx;
                vs[k][j][i] = init((i, j, k), (x, y, z));
                (vs[k][j][i], trs[k][j][i]) = constraints(t, vs[k][j][i]);
            }
        }
    }

    let eigx = Eigenvalues::Analytical(&eigenvaluesx);
    let eigy = Eigenvalues::Analytical(&eigenvaluesy);
    let eigz = Eigenvalues::Analytical(&eigenvaluesz);

    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &ghost, // use noboundary to emulate 1D
        post_constraints: None,
        local_interaction: [1, 1, 1], // use a distance of 0 to emulate 1D
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
        opt: (coord, [eigx, eigy, eigz]),
        p,
        dpde,
        freezeout_energy: None,
    };

    let e = 1e-2;
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_IDEAL_3D]; XY]; XY]; Z],
                   _trs: &[[[[f64; C_IDEAL_3D]; XY]; XY]; Z]| {
        let m = vs
            .iter()
            .flat_map(|v| v.iter().flat_map(|v| v.iter().map(|v| v[0])))
            .sum::<f64>()
            / (XY * XY * Z) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order)
    };

    let observables: [Observable<F_IDEAL_3D, C_IDEAL_3D, XY, XY, Z>; 0] = [];

    run(
        context,
        name,
        crate::hydro::Viscosity::Ideal,
        &names,
        &observables,
        &err_thr,
    )
}
