use crate::{
    hydro::{Viscosity, C_SHEAR_2D, F_SHEAR_2D, HBARC},
    solver::{
        context::{Boundary, Context},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{Eos, Init2D, VOID};

use nalgebra::matrix;

fn gen_constraints<'a>(
    _er: f64,
    p: Eos<'a>,
    dpde: Eos<'a>,
    _temp: Eos<'a>,
) -> Box<dyn Fn(f64, [f64; 6], [f64; 13]) -> ([f64; 6], [f64; 13]) + 'a + Sync> {
    Box::new(
        move |t, [tt00, tt01, tt02, utpi11, utpi12, utpi22], _otrs| {
            // let oe = otrs[0].max(VOID);
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;
            let m = (t01 * t01 + t02 * t02).sqrt();
            let t00 = t00.max(m * (1.0 + 1e-15));

            let sv = |v: f64| m / (t00 + p((t00 - m * v).max(VOID)));
            let eps = 1e-10;
            let mut v = newton(eps, 0.5, |v| sv(v) - v, |v| v.max(0.0).min(1.0));

            // set maximum allowed gamma
            let gamma_max = 20.0f64;
            if (1.0 - v * v) * gamma_max.powi(2) < 1.0 {
                v = (1.0 - 1.0 / gamma_max.powi(2)).sqrt();
                // let e = (t00 - m * v).max(VOID).min(1e10);
                // let t = temp(e) * HBARC;
                // if t > 20.0 {
                //     println!("T: {:e}", t);
                // }
            }

            let e = (t00 - m * v).max(VOID).min(1e10);
            let pe = p(e);
            let vs2 = dpde(e);
            let g = (1.0 - v * v).sqrt();
            // let ut = ((t00 + pe) / (e + pe)).sqrt().max(1.0);
            let ux = g * t01 / (e + pe);
            let uy = g * t02 / (e + pe);
            let ut = (1.0 + ux * ux + uy * uy).sqrt();

            let pi11 = utpi11 / ut;
            let pi12 = utpi12 / ut;
            let pi22 = utpi22 / ut;
            let pi01 = (ux * pi11 + uy * pi12) / ut;
            let pi02 = (ux * pi12 + uy * pi22) / ut;
            let pi00 = (ux * pi01 + uy * pi02) / ut;
            let pi33 = (pi00 - pi01 - pi02) / (t * t); // TODO check division by t*t

            let m = matrix![
                pi11, pi12, 0.0;
                pi12, pi22, 0.0;
                0.0, 0.0, pi33;
            ];
            let eigs = m.symmetric_eigenvalues();
            let smallest = eigs.into_iter().map(|v| *v).min_by(f64::total_cmp).unwrap();

            let epe = e + pe;
            // check that viscosity does not make pressure negative
            let mut r = if epe + smallest < 0.0 {
                -epe / smallest * 0.999
            } else {
                1.0
            };

            while (t00 + r * pi00).powi(2) < (t01 + r * pi01).powi(2) + (t02 + r * pi02).powi(2) {
                r *= 0.99;
            }

            let vs = [
                t00 * t,
                t01 * t,
                t02 * t,
                r * pi11 * ut,
                r * pi12 * ut,
                r * pi22 * ut,
            ];

            let trans = [
                e,
                pe,
                vs2,
                ut,
                ux,
                uy,
                r * pi00,
                r * pi01,
                r * pi02,
                r * pi11,
                r * pi12,
                r * pi22,
                r * pi33,
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
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, pi11, pi12, pi22, _pi33]: [f64; 13],
) -> [f64; 6] {
    [
        t * ((e + pe) * ut * ut - pe), // no pi in Ttt because it is substracted in the flux
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ut * uy),
        ut * pi11,
        ut * pi12,
        ut * pi22,
    ]
}
fn fxuxpi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, pi01, _pi02, pi11, pi12, pi22, _pi33]: [f64; 13],
) -> [f64; 6] {
    [
        t * ((e + pe) * ux * ut + pi01),
        t * ((e + pe) * ux * ux + pe + pi11),
        t * ((e + pe) * ux * uy + pi12),
        ux * pi11,
        ux * pi12,
        ux * pi22,
    ]
}
fn fyuypi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, pi02, pi11, pi12, pi22, _pi33]: [f64; 13],
) -> [f64; 6] {
    [
        t * ((e + pe) * uy * ut + pi02),
        t * ((e + pe) * uy * ux + pi12),
        t * ((e + pe) * uy * uy + pe + pi22),
        uy * pi11,
        uy * pi12,
        uy * pi22,
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
    [_ov, vs]: [&[[[f64; F_SHEAR_2D]; V]; V]; 2],
    [otrs, trs]: [&[[[f64; C_SHEAR_2D]; V]; V]; 2],
    constraints: Constraint<F_SHEAR_2D, C_SHEAR_2D>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    &((etas_min, etas_slope, etas_crv), temperature, tempcut): &((f64, f64, f64), Eos, f64),
) -> [f64; F_SHEAR_2D] {
    let theta = 1.1;

    let pre = &|_t: f64, mut vs: [f64; 6]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        vs[0] = m;
        vs
    };
    let post = &|_t: f64, mut vs: [f64; 6]| {
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
        &u,
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
        &u,
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
    let mev = temp * HBARC;
    let s = (e + pe) / temp;
    let tc = 154.0; // MeV
    let etaovers = (etas_min + etas_slope * (mev - tc) * (mev / tc).powf(etas_crv)).max(etas_min);
    let mut eta = etaovers * s;
    let taupi = 5.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0

    if mev < tempcut {
        eta = 0.0;
    }

    let mut spi = [0.0f64; 6];
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
    }

    let s00 = pe + t * t * pi33;
    [
        -dxf[0] - dyf[0] - dtpi[0] - s00,
        -dxf[1] - dyf[1] - dtpi[1],
        -dxf[2] - dyf[2] - dtpi[2],
        -dxf[3] - dyf[3] + spi[3],
        -dxf[4] - dyf[4] + spi[4],
        -dxf[5] - dyf[5] + spi[5],
    ]
}

pub fn momentum_anysotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &[[[f64; F_SHEAR_2D]; VX]; VY],
    tran: &[[[f64; C_SHEAR_2D]; VX]; VY],
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
    init: Init2D<F_SHEAR_2D>,
    etaovers: (f64, f64, f64),
    shear_temp_cut: f64,
    freezeout_temp_mev: f64,
) -> Option<(
    ([[[f64; F_SHEAR_2D]; V]; V], [[[f64; C_SHEAR_2D]; V]; V]),
    f64,
    usize,
    usize,
)> {
    let constraints = gen_constraints(er, &p, &dpde, temperature);
    let mut vs = [[[0.0; 6]; V]; V];
    let mut trs = [[[0.0; 13]; V]; V];
    let names = (
        ["tt00", "tt01", "tt02", "utpi11", "utpi12", "utpi22"],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "pi00", "pi01", "pi02", "pi11", "pi12", "pi22",
            "pi33",
        ],
    );
    let k = [[[[0.0; 6]; V]; V]; S];
    let v2 = ((V - 1) as f64) / 2.0;
    let mut max_e = 0.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[j][i] = init((i, j), (x, y));
            trs[j][i][0] = vs[j][i][0];
            (vs[j][i], trs[j][i]) = constraints(t, vs[j][i], trs[j][i]);
            max_e = trs[j][i][0].max(max_e);
        }
    }
    // let e = max_e;
    // let t = temperature(e);
    // let s = (e + p(e)) / t;
    // println!("T({}): {}, s: {}", max_e, t * HBARC, s);
    // return None;

    let freezeout_temp = freezeout_temp_mev / HBARC;
    let freezeout_energy = newton(
        1e-10,
        freezeout_temp,
        |e| temperature(e) - freezeout_temp,
        |e| e.max(0.0).min(1000.0),
    );

    let context = Context {
        fun: &flux,
        constraints: &constraints,
        boundary: &[&ghost, &ghost], // use noboundary to emulate 1D
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
        er,
        t,
        ot: t - 1.0,
        t0: t,
        tend,
        opt: (etaovers, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: Some(freezeout_energy),
    };

    let observables: [Observable<F_SHEAR_2D, C_SHEAR_2D, V, V>; 1] =
        [("momentum_anysotropy", &momentum_anysotropy::<V, V>)];

    let temp = shear_temp_cut / HBARC;
    let ecut = newton(
        1e-10,
        temp,
        |e| temperature(e) - temp,
        |e| e.max(0.0).min(1e10),
    );

    run(
        context,
        name,
        Viscosity::Shear(etaovers, ecut),
        &names,
        &observables,
    )
}
