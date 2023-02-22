use crate::{
    hydro::eos::wb,
    solver::{
        context::{Boundary, Context, Integration},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{Eos, Init2D, VOID};

use nalgebra::base::Matrix3;

fn gen_constraints<'a>(
    er: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn(f64, [f64; 10]) -> ([f64; 10], [f64; 13]) + 'a + Sync> {
    Box::new(
        move |t, [tt00, tt01, tt02, utpi00, utpi01, utpi02, utpi11, utpi12, utpi22, mut utpi33]| {
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;
            let m = (t01 * t01 + t02 * t02).sqrt();
            let t00 = t00.max(m * (1.0 + 1e-15));

            let sv = |v: f64| m / (t00 + wb::p((t00 - m * v).max(VOID)));
            let v = newton(1e-5, 0.95, |v| sv(v) - v, |v| v.max(0.0).min(1.0));

            let e = (t00 - m * v).max(VOID);
            let pe = p(e);
            let vs2 = dpde(e);
            let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
            let ux = t01 / ((e + pe) * ut);
            let uy = t02 / ((e + pe) * ut);
            let ut = (1.0 + ux * ux + uy * uy).sqrt();

            let mut utpi = [
                [utpi00, utpi01, utpi02],
                [utpi01, utpi11, utpi12],
                [utpi02, utpi12, utpi22],
            ];
            // utpi[0][0] = utpi[1][1] + utpi[2][2] + utpi33; // traceless

            let m = Matrix3::from_fn(|j, i| utpi[j][i]);
            let eigs = m.symmetric_eigenvalues();
            let smallest = [eigs[0], eigs[1], eigs[2], utpi33]
                .into_iter()
                .min_by(f64::total_cmp)
                .unwrap()
                / ut;

            let epe = e + pe;
            if epe + smallest < 0.0 {
                // check that viscosity does not make pressure negative
                let m = -epe / smallest;
                for j in 0..3 {
                    for i in 0..3 {
                        utpi[j][i] *= m;
                    }
                }
                utpi33 *= m;
            }

            let vs = [
                t00 * t,
                t01 * t,
                t02 * t,
                utpi[0][0],
                utpi[0][1],
                utpi[0][2],
                utpi[1][1],
                utpi[1][2],
                utpi[2][2],
                utpi33,
            ];

            let trans = [
                e,
                pe,
                vs2,
                ut,
                ux,
                uy,
                utpi[0][0] / ut,
                utpi[0][1] / ut,
                utpi[0][2] / ut,
                utpi[1][1] / ut,
                utpi[1][2] / ut,
                utpi[2][2] / ut,
                utpi33 / ut,
            ];
            (vs, trans)
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

pub fn fitutpi(
    t: f64,
    [e, pe, _, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33]: [f64; 13],
) -> [f64; 10] {
    [
        t * ((e + pe) * ut * ut - pe), // no pi in Ttt because it is substracted in the flux
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ut * uy),
        ut * pi00,
        ut * pi01,
        ut * pi02,
        ut * pi11,
        ut * pi12,
        ut * pi22,
        ut * pi33,
    ]
}
fn fxuxpi(
    t: f64,
    [e, pe, _, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33]: [f64; 13],
) -> [f64; 10] {
    [
        t * ((e + pe) * ux * ut + pi01),
        t * ((e + pe) * ux * ux + pe + pi11),
        t * ((e + pe) * ux * uy + pi12),
        ux * pi00,
        ux * pi01,
        ux * pi02,
        ux * pi11,
        ux * pi12,
        ux * pi22,
        ux * pi33,
    ]
}
fn fyuypi(
    t: f64,
    [e, pe, _, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33]: [f64; 13],
) -> [f64; 10] {
    [
        t * ((e + pe) * uy * ut + pi02),
        t * ((e + pe) * uy * ux + pi12),
        t * ((e + pe) * uy * uy + pe + pi22),
        uy * pi00,
        uy * pi01,
        uy * pi02,
        uy * pi11,
        uy * pi12,
        uy * pi22,
        uy * pi33,
    ]
}

fn u(
    _t: f64,
    [_e, _pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33]: [f64; 13],
) -> [f64; 3] {
    [ut, ux, uy]
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; 10]; V]; V]; 2],
    [otrs, trs]: [&[[[f64; 13]; V]; V]; 2],
    constraints: Constraint<10, 13>,
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
        (vs, trs),
        bound,
        pos,
        Dir::X,
        t,
        &fxuxpi,
        (&u, [0, 1, 2]),
        constraints,
        Eigenvalues::Analytical(&eigenvaluesx),
        pre,
        post,
        dx,
        theta,
    );
    let (dyf, dyu) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::Y,
        t,
        &fyuypi,
        (&u, [0, 1, 2]),
        constraints,
        Eigenvalues::Analytical(&eigenvaluesy),
        pre,
        post,
        dx,
        theta,
    );

    let y = bound[1](pos[1], V);
    let x = bound[0](pos[0], V);
    let [_e, _pe, _dpde, out, oux, ouy, opi00, opi01, opi02, _opi11, _opi12, _opi22, _opi33] =
        otrs[y][x];
    let [e, pe, _dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33] = trs[y][x];
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
        let g33 = -1.0; // instead of -1.0 / (t * t) because it is the 'tilde' one
        let pi33_ns = eta * 2.0 * g33 * (u[0] / t - dcuc / 3.0);
        spi[6] = -(pi33 - pi33_ns) / taupi - 1.0 / 3.0 * pi33 * dcuc - pi33 * u[0] / t;
    }

    let s00 = pe + pi33;
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
            mt11 += fxuxpi(t, tran[j][i])[1];
            mt22 += fyuypi(t, tran[j][i])[2];
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
) -> Option<(
    ([[[f64; 10]; V]; V], [[[f64; 13]; V]; V]),
    f64,
    usize,
    usize,
)> {
    let constraints = gen_constraints(er, &p, &dpde);
    let schemename = r.name;
    let mut vs = [[[0.0; 10]; V]; V];
    let mut trs = [[[0.0; 13]; V]; V];
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
            (vs[j][i], trs[j][i]) = constraints(t, init((i, j), (x, y)));
        }
    }
    match r.integration {
        Integration::Explicit => k[S - 1] = vs, // prepare k[S-1] so that it can be use as older time for time derivatives (in this case approximate time derivatives to be zero at initial time)
        Integration::FixPoint(_) => {}
    }
    let integration = r.integration;
    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
        post_constraints: Some(&constraints),
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
