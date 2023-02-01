use crate::solver::{
    context::{Boundary, Context, Integration},
    run,
    space::{id_flux_limiter, order, Dir, Eigenvalues, Order},
    time::{newton::newton, schemes::Scheme},
    utils::{ghost, zeros},
    Transform,
};

use super::{Init2D, Pressure, VOID};

fn constraints(_t: f64, mut vs: [f64; 9]) -> [f64; 9] {
    let tt00 = vs[0];
    let tt01 = vs[1];
    let tt02 = vs[2];
    let tm = (tt01 * tt01 + tt02 * tt02).sqrt();
    let tt00 = tt00.max(tm * (1.0 + 1e-15));
    // TODO: add constraints on shear viscosity
    vs[0] = tt00;
    vs
}

fn gen_transform<'a>(
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn(f64, [f64; 9]) -> [f64; 12] + 'a + Sync> {
    Box::new(
        move |t, [tt00, tt01, tt02, utpi00, utpi01, utpi02, utpi11, utpi12, utpi22]| {
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;
            // the 'i' in it00 stends for 'ideal'
            let sv = |v: f64| {
                let g = (1.0 - v * v).sqrt(); // inverse of ut
                let it00 = t00 - g * utpi00;
                let it01 = t01 - g * utpi01;
                let it02 = t02 - g * utpi02;
                let m = (it01 * it01 + it02 * it02).sqrt();
                let e = (it00 - m * v).max(VOID);
                let v = m / (it00 + p(e));

                v
            };
            let v = newton(er, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
            let g = (1.0 - v * v).sqrt();

            let it00 = t00 - g * utpi00;
            let it01 = t01 - g * utpi01;
            let it02 = t02 - g * utpi02;
            let m = (it01 * it01 + it02 * it02).sqrt();

            let e = (it00 - m * v).max(VOID);
            let pe = p(e);
            let ut = ((it00 + pe) / (e + pe)).sqrt().max(1.0);
            let ux = it01 / ((e + pe) * ut);
            let uy = it02 / ((e + pe) * ut);
            let ut = (1.0 + ux * ux + uy * uy).sqrt();
            [
                e,
                pe,
                dpde(e),
                ut,
                ux,
                uy,
                g * utpi00,
                g * utpi01,
                g * utpi02,
                g * utpi11,
                g * utpi12,
                g * utpi22,
            ]
        },
    )
}

fn eigenvaluesx(_t: f64, [_, _, dpde, ut, ux, _, _, _, _, _, _, _]: [f64; 12]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    (a.abs() + b.sqrt()) / d
}
fn eigenvaluesy(
    t: f64,
    [e, pe, dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22]: [f64; 12],
) -> f64 {
    eigenvaluesx(
        t,
        [e, pe, dpde, ut, uy, ux, pi00, pi01, pi02, pi11, pi12, pi22],
    )
}

pub fn f00(t: f64, [e, pe, _, ut, _, _, pi00, _, _, _, _, _]: [f64; 12]) -> f64 {
    t * ((e + pe) * ut * ut - pe + pi00)
}
pub fn f01(t: f64, [e, pe, _, ut, ux, _, _, pi01, _, _, _, _]: [f64; 12]) -> f64 {
    t * ((e + pe) * ut * ux + pi01)
}
pub fn f02(t: f64, [e, pe, _, ut, _, uy, _, _, pi02, _, _, _]: [f64; 12]) -> f64 {
    t * ((e + pe) * uy * ut + pi02)
}

fn f11(t: f64, [e, pe, _, _, ux, _, _, _, _, pi11, _, _]: [f64; 12]) -> f64 {
    t * ((e + pe) * ux * ux + pe + pi11)
}
fn f12(t: f64, [e, pe, _, _, ux, uy, _, _, _, _, pi12, _]: [f64; 12]) -> f64 {
    t * ((e + pe) * ux * uy + pi12)
}

fn f22(t: f64, [e, pe, _, _, _, uy, _, _, _, _, _, pi22]: [f64; 12]) -> f64 {
    t * ((e + pe) * uy * uy + pe + pi22)
}

fn u1pi00(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u1pi01(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u1pi02(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u1pi11(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u1pi12(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u1pi22(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}

fn u2pi00(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u2pi01(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u2pi02(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u2pi11(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u2pi12(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}
fn u2pi22(
    _t: f64,
    [_e, _pe, _, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22]: [f64; 12],
) -> f64 {
    0.0
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 9]; V]; V]; 2],
    constraints: Transform<9, 9>,
    transform: Transform<9, 12>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [_ot, t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    ord: &Order,
) -> [f64; 9] {
    let theta = 1.1;
    // let theta = 2.0;

    // let pre = &|_t: f64, vs: [f64; 9]| {
    //     let t00 = vs[0];
    //     let t01 = vs[1];
    //     let t02 = vs[2];
    //     let k = t01 * t01 + t02 * t02;
    //     let m = (t00 * t00 - k).sqrt();
    //     [m, t01, t02]
    // };
    // let post = &|_t: f64, vs: [f64; 9]| {
    //     let m = vs[0];
    //     let t01 = vs[1];
    //     let t02 = vs[2];
    //     let k = t01 * t01 + t02 * t02;
    //     let t00 = (m * m + k).sqrt();
    //     [t00, t01, t02]
    // };

    let pre = &id_flux_limiter;
    let post = &id_flux_limiter;

    let diff = order(*ord);
    let divf1 = diff(
        vs,
        bound,
        pos,
        Dir::X,
        t,
        [
            &f01, &f11, &f12, &u1pi00, &u1pi01, &u1pi02, &u1pi11, &u1pi12, &u1pi22,
        ],
        constraints,
        transform,
        Eigenvalues::Analytical(&eigenvaluesx),
        pre,
        post,
        dx,
        theta,
    );
    let divf2 = diff(
        vs,
        bound,
        pos,
        Dir::Y,
        t,
        [
            &f02, &f12, &f22, &u2pi00, &u2pi01, &u2pi02, &u2pi11, &u2pi12, &u2pi22,
        ],
        constraints,
        transform,
        Eigenvalues::Analytical(&eigenvaluesy),
        pre,
        post,
        dx,
        theta,
    );

    let [_e, pe, _dpde, _ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22] =
        transform(t, vs[bound[1](pos[1], V)][bound[0](pos[0], V)]);
    let s = pe;
    [
        -divf1[0] - divf2[0] - s,
        -divf1[1] - divf2[1],
        -divf1[2] - divf2[2],
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

// viscous hydro is in Milne coordinates
pub fn viscoushydro2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Pressure,
    dpde: Pressure,
    init: Init2D<9>,
    space_order: Order,
) -> Option<([[[f64; 9]; V]; V], f64, usize, usize)> {
    let schemename = r.name;
    let mut vs = [[[0.0; 9]; V]; V];
    let names = (
        [
            "tt00", "tt01", "tt02", "utpi00", "utpi01", "utpi02", "utpi11", "utpi12", "utpi22",
        ],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "pi00", "pi01", "pi02", "pi11", "pi12", "pi22",
        ],
    );
    let mut k = [[[[0.0; 9]; V]; V]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init((i, j), (x, y));
        }
    }
    match r.integration {
        Integration::Explicit => k[S - 1] = vs, // prepare k[S-1] so that it can be use as older time for time derivatives (in this case approximate time derivatives to be zero at initial time)
        Integration::FixPoint => {}
    }
    let transform = gen_transform(er, &p, &dpde);
    let integration = r.integration;
    let t00cut: f64 = match space_order {
        Order::Order3 | Order::Order2 => VOID,
        Order::Order3Cut(t00cut) | Order::Order2Cut(t00cut) => t00cut,
    };
    let post = |_t: f64, vs: [f64; 9]| {
        if vs[0] < t00cut {
            [VOID, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        } else {
            vs
        }
    };
    let post: Option<Transform<9, 9>> = match space_order {
        Order::Order3 | Order::Order2 => None,
        Order::Order3Cut(_) | Order::Order2Cut(_) => Some(&post),
    };
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        transform: &transform,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        post_constraints: post,
        local_interaction: [1, 1], // use a distance of 0 to emulate 1D
        vs,
        total_diff_vs: zeros(&vs),
        k,
        r,
        dt: 1e10,
        dx,
        maxdt,
        er,
        t,
        t0: t,
        tend,
        opt: space_order,
        p,
        dpde,
    };
    run(context, name, &schemename, integration, &names)
}
