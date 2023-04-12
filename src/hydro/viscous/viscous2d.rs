use crate::{
    hydro::{Viscosity, C_BOTH_2D, F_BOTH_2D, HBARC},
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
) -> Box<
    dyn Fn(f64, [f64; F_BOTH_2D], [f64; C_BOTH_2D]) -> ([f64; F_BOTH_2D], [f64; C_BOTH_2D])
        + 'a
        + Sync,
> {
    Box::new(
        move |t, cur @ [tt00, tt01, tt02, utpi11, utpi12, utpi22, utppi], _otrs| {
            // let oe = otrs[0].max(VOID);
            let t00 = tt00 / t;
            let mut t01 = tt01 / t;
            let mut t02 = tt02 / t;
            let m = (t01 * t01 + t02 * t02).sqrt();
            if t00 < m {
                let r = t00 / m * (1.0 - 1e-10);
                t01 *= r;
                t02 *= r;
            }
            let t00 = t00.max(m * (1.0 + 1e-15));

            // let sv = |v: f64| {
            //     let g = (1.0 - v * v).sqrt();
            //     let pi00 = g * utpi00;
            //     let pi01 = g * utpi01;
            //     let pi02 = g * utpi02;
            //     let t00 = t00 - pi00;
            //     let t01 = t01 - pi01;
            //     let t02 = t02 - pi02;
            //     let m = (t01 * t01 + t02 * t02).sqrt();
            //     let t00 = t00.max(m * (1.0 + 1e-15));
            //     m / (t00 + p((t00 - m * v).max(VOID)) + (1.0 - v * v).sqrt() * utppi)
            // };
            let sv = |v: f64| m / (t00 + p((t00 - m * v).max(VOID)));
            let cv = |v: f64| v.max(0.0).min(1.0);
            let v = newton(1e-10, 0.5, |v| sv(v) - v, cv);

            let e = (t00 - m * v).max(VOID).min(1e7);
            let g = (1.0 - v * v).sqrt();

            let pe = p(e);
            let mut ux = g * t01 / (e + pe);
            let mut uy = g * t02 / (e + pe);
            let gamma_max2 = 20.0f64.powi(2);
            // let v_max = (1.0 - 1.0 / gamma_max.powi(2)).sqrt();
            let uu = ux * ux + uy * uy;
            if 1.0 + uu > gamma_max2 {
                let r = ((gamma_max2 - 1.0) / uu).sqrt() * 0.999;
                ux *= r;
                uy *= r;
            }
            let ut = (1.0 + ux * ux + uy * uy).sqrt();

            // check that bulk viscosity does not make pressure negative
            let ppi = g * utppi;
            let epe = e + pe;
            let rp = if epe + ppi < 0.0 {
                -epe / ppi * (1.0 - 1e-10)
            } else {
                1.0
            };
            let ppi = rp * ppi;

            let pi11 = utpi11 / ut;
            let pi12 = utpi12 / ut;
            let pi22 = utpi22 / ut;
            let pi01 = (ux * pi11 + uy * pi12) / ut;
            let pi02 = (ux * pi12 + uy * pi22) / ut;
            let pi00 = (ux * pi01 + uy * pi02) / ut;
            let pi33 = (pi00 - pi01 - pi02) / (t * t);

            let m = matrix![
                pi11, pi12, 0.0;
                pi12, pi22, 0.0;
                0.0, 0.0, pi33;
            ];
            let eigs = m.symmetric_eigenvalues();
            let smallest = eigs.into_iter().map(|v| *v).min_by(f64::total_cmp).unwrap();

            // check that shear viscosity does not make pressure negative
            let epeppi = epe + ppi;
            let r = if epeppi + smallest < 0.0 {
                -epeppi / smallest * (1.0 - 1e-10)
            } else {
                1.0
            };
            let pi11 = r * pi11;
            let pi12 = r * pi12;
            let pi22 = r * pi22;
            let pi01 = r * pi01;
            let pi02 = r * pi02;
            let pi00 = r * pi00;
            let pi33 = r * pi33;

            let vs = [
                t00 * t,
                t01 * t,
                t02 * t,
                pi11 * ut,
                pi12 * ut,
                pi22 * ut,
                ppi * ut,
            ];

            let trans = [
                e,
                pe,
                dpde(e),
                ut,
                ux,
                uy,
                pi00,
                pi01,
                pi02,
                pi11,
                pi12,
                pi22,
                pi33,
                ppi,
            ];

            if vs.iter().any(|v| v.is_nan()) || trans.iter().any(|v| v.is_nan()) {
                panic!(
                    "\n\nNaN in constraint\n{:?}\n{:?}\n{} {}\n{:?}\n{:?}\n\n",
                    _otrs, cur, g, v, vs, trans
                );
            }

            (vs, trans)
        },
    )
}

fn eigenvaluesx(_t: f64, [_, _, dpde, ut, ux, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_2D]) -> f64 {
    let vs2 = dpde;
    let a = ut * ux * (1.0 - vs2);
    let b = (ut * ut - ux * ux - (ut * ut - ux * ux - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    let eig = (a.abs() + b.sqrt()) / d;
    eig
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
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, pi11, pi12, pi22, _pi33, ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    [
        t * ((e + pe) * ut * ut - pe), // no pi in Ttt because it is substracted in the flux
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ut * uy),
        ut * pi11,
        ut * pi12,
        ut * pi22,
        ut * ppi,
    ]
}
fn fxuxpi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, pi01, _pi02, pi11, pi12, pi22, _pi33, ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * ux * ut + pi01),
        t * ((e + pe) * ux * ux + pe + pi11),
        t * ((e + pe) * ux * uy + pi12),
        ux * pi11,
        ux * pi12,
        ux * pi22,
        ux * ppi,
    ]
}
fn fyuypi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, pi02, pi11, pi12, pi22, _pi33, ppi]: [f64; C_BOTH_2D],
) -> [f64; F_BOTH_2D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * uy * ut + pi02),
        t * ((e + pe) * uy * ux + pi12),
        t * ((e + pe) * uy * uy + pe + pi22),
        uy * pi11,
        uy * pi12,
        uy * pi22,
        uy * ppi,
    ]
}

fn u(
    _t: f64,
    [_e, _pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi]: [f64;
        C_BOTH_2D],
) -> [f64; 3] {
    [ut, ux, uy]
}

pub enum Coordinate {
    Cartesian,
    Milne,
}

fn flux<const V: usize>(
    [_ov, vs]: [&[[[f64; F_BOTH_2D]; V]; V]; 2],
    [otrs, trs]: [&[[[f64; C_BOTH_2D]; V]; V]; 2],
    constraints: Constraint<F_BOTH_2D, C_BOTH_2D>,
    bound: &[Boundary; 2],
    pos: [i32; 2],
    dx: f64,
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    &((etas_min, etas_slope, etas_crv), (zetas_max, zetas_width, zetas_peak), temperature, tempcut): &(
        (f64, f64, f64),
        (f64, f64, f64),
        Eos,
        f64,
    ),
) -> [f64; F_BOTH_2D] {
    let theta = 1.1;

    let pre = &|_t: f64, mut vs: [f64; F_BOTH_2D]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let k = t01 * t01 + t02 * t02;
        let m = (t00 * t00 - k).sqrt();
        vs[0] = m;
        vs
    };
    let post = &|_t: f64, mut vs: [f64; F_BOTH_2D]| {
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
    let [_e, _pe, _dpde, out, oux, ouy, opi00, opi01, opi02, _opi11, _opi12, _opi22, _opi33, oppi] =
        otrs[y][x];
    let [e, pe, _dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33, ppi] = trs[y][x];
    let pimn = [[pi00, pi01, pi02], [pi01, pi11, pi12], [pi02, pi12, pi22]];

    let g = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    let u = [ut, ux, uy];
    let mut delta = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            delta[a][b] = g[a][b] - u[a] * u[b];
        }
    }
    let ou = [out, oux, ouy];
    let mut odelta = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            odelta[a][b] = g[a][b] - ou[a] * ou[b];
        }
    }

    let dtu = [(ut - out) / cdt, (ux - oux) / cdt, (uy - ouy) / cdt];
    let dtpi = [
        (t * (pi00 - ppi * delta[0][0]) - ot * (opi00 - oppi * odelta[0][0])) / cdt,
        (t * (pi01 - ppi * delta[0][1]) - ot * (opi01 - oppi * odelta[0][1])) / cdt,
        (t * (pi02 - ppi * delta[0][2]) - ot * (opi02 - oppi * odelta[0][2])) / cdt,
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
    let gev = temp * HBARC;
    let s = (e + pe) / temp;
    let tc = 0.154; // GeV

    let vmev = gev.max(tc); // viscous temperature blocked at tc
    let etaovers = etas_min + etas_slope * (vmev - tc) * (vmev / tc).powf(etas_crv);
    let mut eta = etaovers * s;
    let taupi = 5.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0

    let zetaovers = (zetas_max) / (1.0 + ((vmev - zetas_peak) / zetas_width).powi(2));
    let mut zeta = zetaovers * s;
    let tauppi = taupi; // use shear relaxation time for bulk

    // zeta = 0.0;
    if gev < tempcut {
        eta = 0.0;
        zeta = 0.0;
    }

    let mut spi = [0.0f64; 7];
    {
        let mut i = 0;
        // pi^munu
        for a in 0..3 {
            for b in a..3 {
                let pi_ns = eta
                    * ((0..3)
                        .map(|i| delta[a][i] * du[i][b] + delta[b][i] * du[i][a])
                        .sum::<f64>()
                        - 2.0 / 3.0 * delta[a][b] * dcuc);
                let ipi = -1.0 / 3.0 * pimn[a][b] * dcuc
                    - pimn[a][b] * u[0] / t
                    - (0..3)
                        .map(|i| (u[a] * pimn[b][i] + u[b] * pimn[a][i]) * ddu[i] * g[i][i])
                        .sum::<f64>();
                spi[i] = -(pimn[a][b] - pi_ns) / taupi + ipi;
                i += 1;
            }
        }
        // Pi
        let ppi_ns = -zeta * dcuc;
        let ippi = -1.0 / 3.0 * ppi * dcuc - ppi * u[0] / t;
        spi[i] = -(ppi - ppi_ns) / tauppi + ippi;
    }

    let s00 = (pe + ppi) + t * t * pi33;
    [
        -dxf[0] - dyf[0] - dtpi[0] - s00,
        -dxf[1] - dyf[1] - dtpi[1],
        -dxf[2] - dyf[2] - dtpi[2],
        -dxf[3] - dyf[3] + spi[3],
        -dxf[4] - dyf[4] + spi[4],
        -dxf[5] - dyf[5] + spi[5],
        -dxf[6] - dyf[6] + spi[6],
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
pub fn viscous2d<const V: usize, const S: usize>(
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
    let constraints = gen_constraints(er, &p, &dpde, temperature);

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
            (vs[j][i], trs[j][i]) = constraints(t, vs[j][i], trs[j][i]);
            max_e = trs[j][i][0].max(max_e);
        }
    }
    // let e = max_e;
    // let t = temperature(e);
    // let s = (e + p(e)) / t;
    // println!("T({}): {}, s: {}", max_e, t * HBARC, s);
    // return None;

    let freezeout_temp_fm = freezeout_temp / HBARC;
    let freezeout_energy = newton(
        1e-10,
        freezeout_temp_fm,
        |e| temperature(e) - freezeout_temp_fm,
        |e| e.max(0.0).min(1e10),
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
        opt: (etaovers, zetaovers, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: Some(freezeout_energy),
    };

    let observables: [Observable<F_BOTH_2D, C_BOTH_2D, V, V>; 1] =
        [("momentum_anysotropy", &momentum_anysotropy::<V, V>)];

    let temp_fm = shear_temp_cut / HBARC;
    let ecut = newton(
        1e-10,
        temp_fm,
        |e| temperature(e) - temp_fm,
        |e| e.max(0.0).min(1e10),
    );

    run(
        context,
        name,
        Viscosity::Both(etaovers, zetaovers, ecut),
        &names,
        &observables,
    )
}
