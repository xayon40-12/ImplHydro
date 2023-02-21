use crate::{
    hydro::eos::wb,
    solver::{
        context::{Boundary, Context, Integration},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Observable, Transform,
    },
};

use crate::hydro::{Eos, Init2D, VOID};

use nalgebra::base::Matrix3;

fn constraints(t: f64, mut vs: [f64; 10]) -> [f64; 10] {
    let tt00 = vs[0];
    let tt01 = vs[1];
    let tt02 = vs[2];
    let tm = (tt01 * tt01 + tt02 * tt02).sqrt();
    let tt00 = tt00.max(tm * (1.0 + 1e-15));
    vs[0] = tt00;
    let utpi00 = vs[3];
    let utpi01 = vs[4];
    let utpi02 = vs[5];
    let utpi11 = vs[6];
    let utpi12 = vs[7];
    let utpi22 = vs[8];
    let utpi33 = vs[9];
    vs[3] = utpi11 + utpi22 + t * t * utpi33; // traceless

    const ORTH: bool = false;
    const PRESSURE: bool = true;
    if ORTH || PRESSURE {
        let t00 = tt00 / t;
        let t01 = tt01 / t;
        let t02 = tt02 / t;

        let sv = |v: f64| {
            let m = (t01 * t01 + t02 * t02).sqrt();
            let e = (t00 - m * v).max(VOID);
            let v = m / (t00 + wb::p(e));

            v
        };
        let v = newton(1e-5, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
        let m = (t01 * t01 + t02 * t02).sqrt();

        let e = (t00 - m * v).max(VOID);
        let pe = wb::p(e);
        let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
        let ux = t01 / ((e + pe) * ut);
        let uy = t02 / ((e + pe) * ut);
        let ut = (1.0 + ux * ux + uy * uy).sqrt();

        let pi00 = utpi00 / ut;
        let pi01 = utpi01 / ut;
        let pi02 = utpi02 / ut;
        let pi11 = utpi11 / ut;
        let pi12 = utpi12 / ut;
        let pi22 = utpi22 / ut;
        let mut pi33 = utpi33 / ut;
        let mut pi = [[pi00, pi01, pi02], [pi01, pi11, pi12], [pi02, pi12, pi22]];

        let m = Matrix3::from_fn(|j, i| pi[j][i]);
        let eigs = m.symmetric_eigenvalues();
        let smallest = [eigs[0], eigs[1], eigs[2], t * t * pi33]
            .into_iter()
            .min_by(f64::total_cmp)
            .unwrap();

        if pe + smallest < 0.0 {
            // check that viscosity does not make pressure negative
            let m = -pe / smallest * 0.95;
            for j in 0..3 {
                for i in 0..3 {
                    pi[j][i] *= m;
                }
            }
            pi33 *= m;
        }

        vs[3] = pi00 * ut;
        vs[4] = pi01 * ut;
        vs[5] = pi02 * ut;
        vs[6] = pi11 * ut;
        vs[7] = pi12 * ut;
        vs[8] = pi22 * ut;
        vs[9] = pi33 * ut;

        if ORTH {
            // let utpi00 = ux * pi01 + uy * pi02;
            let utpi01 = ux * pi11 + uy * pi12;
            let utpi02 = ux * pi12 + uy * pi22;
            // if (vs[3] - utpi00).abs() > 1e-6 {
            //     panic!("\nutpi00\n{:e} {:e}", vs[3], utpi00);
            // }
            // vs[3] = utpi00;
            vs[4] = utpi01;
            vs[5] = utpi02;
        }
    }

    vs
}

fn gen_transform<'a>(
    er: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn(f64, [f64; 10]) -> [f64; 13] + 'a + Sync> {
    Box::new(
        move |t, [tt00, tt01, tt02, utpi00, utpi01, utpi02, utpi11, utpi12, utpi22, utpi33]| {
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;

            let sv = |v: f64| {
                let m = (t01 * t01 + t02 * t02).sqrt();
                let e = (t00 - m * v).max(VOID);
                let v = m / (t00 + p(e));

                v
            };
            let v = newton(er, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));
            let m = (t01 * t01 + t02 * t02).sqrt();

            let e = (t00 - m * v).max(VOID);
            let pe = p(e);
            let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
            let ux = t01 / ((e + pe) * ut);
            let uy = t02 / ((e + pe) * ut);
            let ut = (1.0 + ux * ux + uy * uy).sqrt();
            let vs2 = dpde(e);
            [
                e,
                pe,
                vs2,
                ut,
                ux,
                uy,
                utpi00 / ut,
                utpi01 / ut,
                utpi02 / ut,
                utpi11 / ut,
                utpi12 / ut,
                utpi22 / ut,
                utpi33 / ut,
            ]
        },
    )
}

fn eigenvaluesx(_t: f64, [_, _, dpde, ut, ux, _, _, _, _, _, _, _, _]: [f64; 13]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    let eig = (a.abs() + b.sqrt()) / d;
    eig
}
fn eigenvaluesy(
    t: f64,
    [e, pe, dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33]: [f64; 13],
) -> f64 {
    eigenvaluesx(
        t,
        [
            e, pe, dpde, ut, uy, ux, pi00, pi01, pi02, pi11, pi12, pi22, pi33,
        ],
    )
}

pub fn fi00(t: f64, [e, pe, _, ut, _, _, _, _, _, _, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * ut * ut - pe)
}
pub fn fi01(t: f64, [e, pe, _, ut, ux, _, _, _, _, _, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * ut * ux)
}
pub fn fi02(t: f64, [e, pe, _, ut, _, uy, _, _, _, _, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * uy * ut)
}

fn f01(t: f64, [e, pe, _, ut, ux, _, _, pi01, _, _, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * ut * ux + pi01)
}
fn f02(t: f64, [e, pe, _, ut, _, uy, _, _, pi02, _, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * uy * ut + pi02)
}

fn f11(t: f64, [e, pe, _, _, ux, _, _, _, _, pi11, _, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * ux * ux + pe + pi11)
}
fn f12(t: f64, [e, pe, _, _, ux, uy, _, _, _, _, pi12, _, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * ux * uy + pi12)
}

fn f22(t: f64, [e, pe, _, _, _, uy, _, _, _, _, _, pi22, _]: [f64; 13]) -> f64 {
    t * ((e + pe) * uy * uy + pe + pi22)
}

pub fn u0pi00(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi00
}
pub fn u0pi01(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi01
}
pub fn u0pi02(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi02
}
pub fn u0pi11(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, _pi02, pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi11
}
pub fn u0pi12(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi12
}
pub fn u0pi22(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, pi22, _pi33]: [f64; 13],
) -> f64 {
    ut * pi22
}
pub fn u0pi33(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, pi33]: [f64; 13],
) -> f64 {
    ut * pi33
}

fn u1pi00(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi00
}
fn u1pi01(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi01
}
fn u1pi02(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi02
}
fn u1pi11(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, _pi02, pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi11
}
fn u1pi12(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, _pi02, _pi11, pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi12
}
fn u1pi22(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, pi22, _pi33]: [f64; 13],
) -> f64 {
    ux * pi22
}
fn u1pi33(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, pi33]: [f64; 13],
) -> f64 {
    ux * pi33
}

fn u2pi00(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi00
}
fn u2pi01(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi01
}
fn u2pi02(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi02
}
fn u2pi11(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, _pi02, pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi11
}
fn u2pi12(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, _pi02, _pi11, pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi12
}
fn u2pi22(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, pi22, _pi33]: [f64; 13],
) -> f64 {
    uy * pi22
}
fn u2pi33(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, pi33]: [f64; 13],
) -> f64 {
    uy * pi33
}

fn u0(
    _t: f64,
    [_e, _pe, _, ut, _ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ut
}
fn u1(
    _t: f64,
    [_e, _pe, _, _ut, ux, _uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    ux
}
fn u2(
    _t: f64,
    [_e, _pe, _, _ut, _ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> f64 {
    uy
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [ov, vs]: [&[[[f64; 10]; V]; V]; 2],
    constraints: Transform<10, 10>,
    transform: Transform<10, 13>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    (etaovers, temperature): &(f64, Eos),
) -> [f64; 10] {
    let theta = 1.1;

    let pre = &|_t: f64, mut vs: [f64; 10]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        vs[0] = m;
        vs
    };
    let post = &|_t: f64, mut vs: [f64; 10]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let t00 = (m * m + k).sqrt();
        vs[0] = t00;
        vs
    };

    // let pre = &id_flux_limiter;
    // let post = &id_flux_limiter;

    let diff = kt;
    let (dxf, dxu) = diff(
        vs,
        bound,
        pos,
        Dir::X,
        t,
        [
            &f01, &f11, &f12, &u1pi00, &u1pi01, &u1pi02, &u1pi11, &u1pi12, &u1pi22, &u1pi33,
        ],
        ([&u0, &u1, &u2], [0, 1, 2]),
        constraints,
        transform,
        Eigenvalues::Analytical(&eigenvaluesx),
        pre,
        post,
        dx,
        theta,
    );
    let (dyf, dyu) = diff(
        vs,
        bound,
        pos,
        Dir::Y,
        t,
        [
            &f02, &f12, &f22, &u2pi00, &u2pi01, &u2pi02, &u2pi11, &u2pi12, &u2pi22, &u2pi33,
        ],
        ([&u0, &u1, &u2], [0, 1, 2]),
        constraints,
        transform,
        Eigenvalues::Analytical(&eigenvaluesy),
        pre,
        post,
        dx,
        theta,
    );

    let [_e, _pe, _dpde, out, oux, ouy, opi00, opi01, opi02, _opi11, _opi12, _opi22, _opi33] =
        transform(
            t,
            constraints(t, ov[bound[1](pos[1], V)][bound[0](pos[0], V)]),
        );
    let [e, pe, _dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33] = transform(
        t,
        constraints(t, vs[bound[1](pos[1], V)][bound[0](pos[0], V)]),
    );
    let pi = [[pi00, pi01, pi02], [pi01, pi11, pi12], [pi02, pi12, pi22]];

    let u = [ut, ux, uy];
    let g = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    let mut delta = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            delta[a][b] = g[a][b] - u[a] * u[b];
        }
    }

    let dtu = [(ut - out) / cdt, (ux - oux) / cdt, (uy - ouy) / cdt];
    let dtpi = [
        (t * pi00 - ot * opi00) / cdt,
        (t * pi01 - ot * opi01) / cdt,
        (t * pi02 - ot * opi02) / cdt,
    ];
    let du = [dtu, dxu, dyu];
    let mut ddu = [0.0f64; 3];
    let mut dcuc = u[0] / t;
    for j in 0..3 {
        dcuc += du[j][j];
        for i in 0..3 {
            ddu[j] += u[i] * du[i][j];
        }
    }

    let temp = temperature(e);
    let mev = temp * 197.3;
    let cutmev = 50.0;
    let s = (e + pe) / temp;
    let mut eta = etaovers * s;
    let taupi = 3.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0
    if mev < cutmev {
        eta = 0.0;
    }

    let mut spi = [0.0f64; 7];
    {
        let mut i = 0;
        for a in 0..3 {
            for b in a..3 {
                let pi_ns: f64 = eta
                    * ((0..3)
                        .map(|i| delta[a][i] * du[i][b] + delta[b][i] * du[i][a])
                        .sum::<f64>()
                        - 2.0 / 3.0 * delta[a][b] * dcuc);
                let ipi = -1.0 / 3.0 * pi[a][b] * dcuc
                    - pi[a][b] * u[0] / t
                    - (0..3)
                        .map(|i| (u[a] * pi[b][i] + u[b] * pi[a][i]) * ddu[i] * g[i][i])
                        .sum::<f64>();
                spi[i] = -(pi[a][b] - pi_ns) / taupi + ipi;
                i += 1;
            }
        }
        let g33 = -1.0 / (t * t);
        let pi33_ns = eta * 2.0 * g33 * (u[0] / t - dcuc / 3.0);
        spi[6] = -(pi33 - pi33_ns) / taupi - 1.0 / 3.0 * pi33 * dcuc - pi33 * u[0] / t;
    }

    let s00 = pe + t * t * pi33;
    [
        -dxf[0] - dyf[0] - dtpi[0] - s00,
        -dxf[1] - dyf[1] - dtpi[1],
        -dxf[2] - dyf[2] - dtpi[2],
        -dxf[3] - dyf[3] + spi[0],
        -dxf[4] - dyf[4] + spi[1],
        -dxf[5] - dyf[5] + spi[2],
        -dxf[6] - dyf[6] + spi[3],
        -dxf[7] - dyf[7] + spi[4],
        -dxf[8] - dyf[8] + spi[5],
        -dxf[9] - dyf[9] + spi[6],
    ]
}

pub fn momentum_anysotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &[[[f64; 10]; VX]; VY],
    tran: &[[[f64; 13]; VX]; VY],
) -> Vec<f64> {
    let mut mt11 = 0.0;
    let mut mt22 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            mt11 += f11(t, tran[j][i]);
            mt22 += f22(t, tran[j][i]);
        }
    }
    let anysotropy = (mt11 - mt22) / (mt11 + mt22);
    vec![anysotropy]
}

// viscous hydro is in Milne coordinates
pub fn shear2d<const V: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    er: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    temperature: Eos,
    init: Init2D<10>,
    etaovers: f64,
) -> Option<([[[f64; 10]; V]; V], f64, usize, usize)> {
    let schemename = r.name;
    let mut vs = [[[0.0; 10]; V]; V];
    let names = (
        [
            "tt00", "tt01", "tt02", "utpi00", "utpi01", "utpi02", "utpi11", "utpi12", "utpi22",
            "utpi33",
        ],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "pi00", "pi01", "pi02", "pi11", "pi12", "pi22",
            "pi33",
        ],
    );
    let mut k = [[[[0.0; 10]; V]; V]; S];
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
        Integration::FixPoint(_) => {}
    }
    let transform = gen_transform(er, &p, &dpde);
    let integration = r.integration;
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        transform: &transform,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        post_constraints: Some(&constraints),
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
        opt: (etaovers, temperature),
        p,
        dpde,
    };

    let observables: [Observable<10, 13, V, V>; 1] =
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
