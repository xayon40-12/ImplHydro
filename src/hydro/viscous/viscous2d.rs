use crate::{
    hydro::{utils::eigenvaluesk, Viscosity, C_MILNE_BOTH_2D, FREESTREAM_2D, F_BOTH_2D, HBARC},
    solver::{
        context::{Arr, BArr, Boundary, Context, Integration, DIM},
        run,
        space::{kt::kt, Eigenvalues, FluxInfo, InDir::*},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable, EXACT,
    },
};
use boxarray::boxarray;

use crate::hydro::{Eos, Init2D, VOID};

use nalgebra::matrix;

pub fn init_from_entropy_density_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    s: &'a [[f64; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
    _entropy: Eos<'a>,
    temperature: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_BOTH_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let s = s[j][i] / HBARC / t0;
        // let e = newton(1e-10, s, |e| entropy(e) - s, |e| e.max(0.0).min(1e10)).max(VOID);
        let e = s;
        let vars = [
            e,
            p(e),
            dpde(e),
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            temperature(e),
        ];
        fitutpi(t0, vars)
    })
}

pub fn init_from_freestream_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    trs: &'a [[[f64; FREESTREAM_2D]; VX]; VY],
    p: Eos<'a>,
    dpde: Eos<'a>,
    temperature: Eos<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; F_BOTH_2D] + 'a> {
    Box::new(move |(i, j), _| {
        let e = trs[j][i][0].max(VOID);
        let ut = trs[j][i][1];
        let ux = trs[j][i][2];
        let uy = trs[j][i][3];
        // TODO use viscosity from freestream
        // let pi00 = trs[j][i][4];
        // let pi01 = trs[j][i][5];
        // let pi02 = trs[j][i][6];
        // let pi11 = trs[j][i][7];
        // let pi12 = trs[j][i][8];
        // let pi22 = trs[j][i][9];
        // let pi33 = pi00 - pi11 - pi22;
        // let bulk = trs[j][i][10];
        let vars = [
            e,
            p(e),
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
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            temperature(e),
            // pi00,
            // pi01,
            // pi02,
            // pi11,
            // pi12,
            // pi22,
            // pi33,
            // bulk,
        ];
        fitutpi(t0, vars)
    })
}

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    temp: Eos<'a>,
    _implicit: bool,
) -> Box<dyn Fn(f64, [f64; F_BOTH_2D]) -> ([f64; F_BOTH_2D], [f64; C_MILNE_BOTH_2D]) + 'a + Sync> {
    Box::new(
        move |t, [tt00, tt01, tt02, utpi11, utpi12, utpi22, utppi]| {
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;
            let m = (t01 * t01 + t02 * t02).sqrt();
            let t00 = t00.max(m * (1.0 + 1e-15));

            let sv = |v: f64| m / (t00 + p((t00 - m * v).max(VOID)));
            let cv = |v: f64| v.max(0.0).min(1.0);
            let v = newton(1e-10, 0.5, |v| sv(v) - v, cv);

            let e = (t00 - m * v).max(VOID).min(1e10);

            {
                let g = (1.0 - v * v).sqrt();

                let pe = p(e);
                let ux = g * t01 / (e + pe);
                let uy = g * t02 / (e + pe);
                let ut = (1.0 + ux * ux + uy * uy).sqrt();

                // check that bulk viscosity does not make pressure negative
                let ppi = g * utppi;
                let epe = e + pe;

                let pi11 = utpi11 / ut;
                let pi12 = utpi12 / ut;
                let pi22 = utpi22 / ut;
                let pi01 = (ux * pi11 + uy * pi12) / ut;
                let pi02 = (ux * pi12 + uy * pi22) / ut;
                let pi00 = (ux * pi01 + uy * pi02) / ut;
                let pi33 = pi00 - pi11 - pi22;

                // let m = matrix![ // \pi^{\mu\nu}   \pi^\mu_\nu
                //     pi00, pi01, pi02;
                //     pi01, pi11, pi12;
                //     pi02, pi12, pi22;
                // ]; // FIXME this is $\pi^{\mu\nu}$ but we want the eigenvalues of $\pi^\mu_\nu$
                let m = matrix![ // \pi^{\mu\nu}   \pi^\mu_\nu
                    pi11, pi12;
                    pi12, pi22;
                ]; // FIXME this is $\pi^{\mu\nu}$ but we want the eigenvalues of $\pi^\mu_\nu$
                let mut eigs: Vec<f64> = m
                    .symmetric_eigenvalues()
                    .into_iter()
                    .cloned()
                    .chain([0.0, pi33].into_iter())
                    .collect();
                eigs.sort_by(f64::total_cmp);
                let smallest = eigs[0];

                // check that shear viscosity does not make pressure negative
                // let r = if epe + ppi + smallest < 0.0 {
                let r = if epe + ppi + smallest < 0.0 {
                    -epe / (ppi + smallest) * (1.0 - 1e-10)
                } else {
                    1.0
                };

                let ppi = r * ppi;
                let pi11 = r * pi11;
                let pi12 = r * pi12;
                let pi22 = r * pi22;
                let pi01 = r * pi01;
                let pi02 = r * pi02;
                let pi00 = r * pi00;
                let pi33 = r * pi33;
                let lam0 = r * eigs[0];
                let lam1 = r * eigs[1];
                let lam2 = r * eigs[2];
                let lam3 = r * eigs[3];

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
                    lam0,
                    lam1,
                    lam2,
                    lam3,
                    r,
                    temp(e),
                ];

                // if vs.iter().any(|v| v.is_nan()) || trans.iter().any(|v| v.is_nan()) {
                //     panic!(
                //         "\n\nNaN in constraint\n{:?}\n{} {}\n{:?}\n{:?}\n\n",
                //         cur, g, v, vs, trans
                //     );
                // }

                (vs, trans)
            }
        },
    )
}

fn eigenvaluesx(
    _t: f64,
    [_, _, dpde, ut, ux, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _]: [f64; C_MILNE_BOTH_2D],
) -> f64 {
    eigenvaluesk(dpde, ut, ux)
}
fn eigenvaluesy(
    _t: f64,
    [_e, _pe, dpde, ut, _ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi, _, _, _, _, _, _]: [f64;
        C_MILNE_BOTH_2D],
) -> f64 {
    eigenvaluesk(dpde, ut, uy)
}

pub fn fitutpi(
    t: f64,
    [e, pe, _, ut, ux, uy, _pi00, _pi01, _pi02, pi11, pi12, pi22, _pi33, ppi, _, _, _, _, _, _]: [f64;
        C_MILNE_BOTH_2D],
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
    [e, pe, _, ut, ux, uy, _pi00, pi01, _pi02, pi11, pi12, pi22, _pi33, ppi, _, _, _, _, _, _]: [f64;
        C_MILNE_BOTH_2D],
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
    [e, pe, _, ut, ux, uy, _pi00, _pi01, pi02, pi11, pi12, pi22, _pi33, ppi, _, _, _, _, _, _]: [f64;
        C_MILNE_BOTH_2D],
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
    [_e, _pe, _, ut, ux, uy, _pi00, _pi01, _pi02, _pi11, _pi12, _pi22, _pi33, _ppi, _, _, _, _, _, _]: [f64;
        C_MILNE_BOTH_2D],
) -> [f64; 3] {
    [ut, ux, uy]
}

fn flux<const V: usize>(
    k: &Arr<F_BOTH_2D, V, V, 1>,
    [_ov, vs]: [&Arr<F_BOTH_2D, V, V, 1>; 2],
    [otrs, trs]: [&Arr<C_MILNE_BOTH_2D, V, V, 1>; 2],
    constraints: Constraint<F_BOTH_2D, C_MILNE_BOTH_2D>,
    bound: Boundary<F_BOTH_2D, V, V, 1>,
    pos: [i32; DIM],
    dxs: [f64; DIM],
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    &(
        (etas_min, etas_slope, etas_crv),
        (zetas_max, zetas_width, zetas_peak),
        entropy,
        temperature,
        _tempcut,
    ): &((f64, f64, f64), (f64, f64, f64), Eos, Eos, f64),
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

    let diff = kt;
    let flux_infos = [
        X(FluxInfo {
            flux: &fxuxpi,
            secondary: &u,
            eigenvalues: Eigenvalues::Analytical(&eigenvaluesx),
        }),
        Y(FluxInfo {
            flux: &fyuypi,
            secondary: &u,
            eigenvalues: Eigenvalues::Analytical(&eigenvaluesy),
        }),
    ];
    let [(dxf, dxu), (dyf, dyu)] = diff(
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

    let z = 0;
    let y = pos[1] as usize;
    let x = pos[0] as usize;
    let [ktt00, ktt01, ktt02, kutpi11, kutpi12, kutpi22, kutppi] = k[z][y][x];
    let [tt00, tt01, tt02, _utpi11, _utpi12, _utpi22, _utppi] = vs[z][y][x];
    let [_oe, _ope, _odpde, out, oux, ouy, opi00, opi01, opi02, _opi11, _opi12, _opi22, _opi33, oppi, _, _, _, _, _, _] =
        otrs[z][y][x];
    let [e, pe, dpde, ut, ux, uy, pi00, pi01, pi02, pi11, pi12, pi22, pi33, ppi, _, _, _, _, _, _] =
        trs[z][y][x];
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

    let (dtu, dtpi) = {
        if EXACT {
            let ut2 = ut * ut;
            let ux2 = ux * ux;
            let uy2 = uy * uy;
            let d0 = dpde * ut2 - dpde + ut2;
            let dp0 = dpde * ut2 + dpde + ut2;
            let d1 = d0 * ut2 - dp0 * ux2;
            let d2 = d0 * ut2 - dp0 * uy2;
            let d12 = d0 * ut2 - dp0 * ux2 - dp0 * uy2;
            let ep = e + pe;
            let u0ep = ut * ep;
            let dpde1 = dpde + 1.0;

            // V = (\epsilon, u^x, u^y)
            // v = dV/dt
            // k = d T^{t\nu}/dt = dk/dv * v    because k is linear in d/dt
            // M = dk/dv
            // k = M*v
            // v = M^-1 k
            // m = M^-1
            // v = m k
            //
            //     ⎡    2     2     2                                        ⎤
            //     ⎢  u₀  + u₁  + u₂         -2⋅u₀⋅u₁           -2⋅u₀⋅u₂     ⎥
            //     ⎢  ───────────────        ─────────          ─────────    ⎥
            //     ⎢        D₁₂                 D₁₂                D₁₂       ⎥
            //     ⎢                                                         ⎥
            //     ⎢   2                                                     ⎥
            //     ⎢-u₀ ⋅u₁⋅(dpde + 1)           D₂             DP₀⋅u₁⋅u₂    ⎥
            // m = ⎢───────────────────  ─────────────────  ─────────────────⎥
            //     ⎢   D₁₂⋅(e + p(e))    D₁₂⋅u₀⋅(e + p(e))  D₁₂⋅u₀⋅(e + p(e))⎥
            //     ⎢                                                         ⎥
            //     ⎢   2                                                     ⎥
            //     ⎢-u₀ ⋅u₂⋅(dpde + 1)       DP₀⋅u₁⋅u₂              D₁       ⎥
            //     ⎢───────────────────  ─────────────────  ─────────────────⎥
            //     ⎣   D₁₂⋅(e + p(e))    D₁₂⋅u₀⋅(e + p(e))  D₁₂⋅u₀⋅(e + p(e))⎦
            let d12m = [
                [ut2 + ux2 + uy2, -2.0 * ut * ux, -2.0 * ut * uy],
                [-ut2 * ux * dpde1 / ep, d2 / u0ep, dp0 * ux * uy / u0ep],
                [-ut2 * uy * dpde1 / ep, dp0 * ux * uy / u0ep, d1 / u0ep],
            ];
            // d_t (tT^{t\nu}) = T^{t\nu} + t d_t T^{t\nu}
            let k = [
                (ktt00 - tt00 / t) / t,
                (ktt01 - tt01 / t) / t,
                (ktt02 - tt02 / t) / t,
            ];
            let mut v = [0.0f64; 3];
            for j in 0..3 {
                for i in 0..3 {
                    v[j] += d12m[j][i] / d12 * k[i];
                }
            }
            let _dedt = v[0];
            let dtu = [(ux * v[1] + uy * v[2]) / ut, v[1], v[2]];
            // d_t(ut*Pi) = Pi*d_t ut + ut*d_t Pi
            let dtppi = (kutppi - ppi * dtu[0]) / ut;
            // d_t(ut*pi^{t\nu}) = pi^{t\nu}*d_t ut + ut*d_t pi^{t\nu}
            let dtpi11 = (kutpi11 - pimn[1][1] * dtu[0]) / ut;
            let dtpi12 = (kutpi12 - pimn[1][2] * dtu[0]) / ut;
            let dtpi22 = (kutpi22 - pimn[2][2] * dtu[0]) / ut;
            // utpi10 = uxpi11 + uypi12
            // kutpi10 = uxdtpi11 + pi11dtux + uydtpi12 + pi12dtuy
            let kutpi01 = ux * dtpi11 + pi11 * dtu[1] + uy * dtpi12 + pi12 * dtu[2];
            let dtpi01 = (kutpi01 - pimn[1][0] * dtu[0]) / ut;
            // utpi20 = uxpi21 + uypi22
            let kutpi02 = ux * dtpi12 + pi12 * dtu[1] + uy * dtpi22 + pi22 * dtu[2];
            let dtpi02 = (kutpi02 - pimn[2][0] * dtu[0]) / ut;
            // utpi00 = uxpi10 + uypi20
            let kutpi00 = ux * dtpi01 + pi01 * dtu[1] + uy * dtpi02 + pi02 * dtu[2];
            let dtpi00 = (kutpi00 - pimn[0][0] * dtu[0]) / ut;
            let dttpi00 = pi00 + t * dtpi00;
            let dttpi01 = pi01 + t * dtpi01;
            let dttpi02 = pi02 + t * dtpi02;
            let mut dttppideltat = [0.0f64; 3];
            for i in 0..3 {
                dttppideltat[i] = ppi * delta[0][i]
                    + t * (delta[0][i] * dtppi - ppi * (u[i] * dtu[0] + u[0] * dtu[i]));
            }
            let dtpi = [
                dttpi00 - dttppideltat[0],
                dttpi01 - dttppideltat[1],
                dttpi02 - dttppideltat[2],
            ];

            (dtu, dtpi)
        } else {
            let dtu = [(ut - out) / cdt, (ux - oux) / cdt, (uy - ouy) / cdt];
            let dtpi = [
                (t * (pi00 - ppi * delta[0][0]) - ot * (opi00 - oppi * odelta[0][0])) / cdt,
                (t * (pi01 - ppi * delta[0][1]) - ot * (opi01 - oppi * odelta[0][1])) / cdt,
                (t * (pi02 - ppi * delta[0][2]) - ot * (opi02 - oppi * odelta[0][2])) / cdt,
            ];
            (dtu, dtpi)
        }
    };
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
    let s = entropy(e); // (e + pe) / temp;
    let tc = 0.154; // GeV

    let vgev = gev.max(tc); // viscous temperature [GeV] blocked at tc
    let etaovers = etas_min + etas_slope * (vgev - tc) * (vgev / tc).powf(etas_crv);
    let eta = etaovers * s;
    let taupi = 5.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0

    let zetaovers = (zetas_max) / (1.0 + ((vgev - zetas_peak) / zetas_width).powi(2));
    let zeta = zetaovers * s;
    let tauppi = taupi; // use shear relaxation time for bulk

    // if gev < tempcut {
    //     let tau_decay = 10.0;
    //     let m = ((1.0 - tempcut / gev) / tau_decay).exp();
    //     eta *= m;
    //     zeta *= m;
    // }

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

    let s00 = (pe + ppi) + pi33;
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

pub fn momentum_anisotropy<const VX: usize, const VY: usize>(
    t: f64,
    _vs: &[[[[f64; F_BOTH_2D]; VX]; VY]; 1],
    tran: &[[[[f64; C_MILNE_BOTH_2D]; VX]; VY]; 1],
) -> Vec<f64> {
    let mut mt11 = 0.0;
    let mut mt12 = 0.0;
    let mut mt22 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            let [_, t11, t12, ..] = fxuxpi(t, tran[0][j][i]);
            let [_, _, t22, ..] = fyuypi(t, tran[0][j][i]);
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

// viscous hydro is in Milne coordinates
pub fn viscous2d<const V: usize, const S: usize>(
    name: &(&str, usize),
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    entropy: Eos,
    temperature: Eos,
    init: Init2D<F_BOTH_2D>,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    shear_temp_cut: f64,
    freezeout_temp: f64,
    save_raw: Option<f64>,
) -> Option<(
    (BArr<F_BOTH_2D, V, V, 1>, BArr<C_MILNE_BOTH_2D, V, V, 1>),
    f64,
    usize,
    usize,
)> {
    let implicit = match r.integration {
        Integration::FixPoint => true,
        Integration::Explicit => false,
    };
    let constraints = gen_constraints(&p, &dpde, temperature, implicit);

    let mut vs: Box<[[[[f64; F_BOTH_2D]; V]; V]; 1]> = boxarray(0.0);
    let mut trs: Box<[[[[f64; C_MILNE_BOTH_2D]; V]; V]; 1]> = boxarray(0.0);
    let names = (
        ["tt00", "tt01", "tt02", "utpi11", "utpi12", "utpi22", "utPi"],
        [
            "e",
            "pe",
            "dpde",
            "ut",
            "ux",
            "uy",
            "pi00",
            "pi01",
            "pi02",
            "pi11",
            "pi12",
            "pi22",
            "pi33",
            "Pi",
            "lam0",
            "lam1",
            "lam2",
            "lam3",
            "renorm",
            "temperature",
        ],
    );
    let k: Box<[[[[[f64; F_BOTH_2D]; V]; V]; 1]; S]> = boxarray(0.0);
    let v2 = ((V - 1) as f64) / 2.0;
    let mut max_e = 0.0;
    for j in 0..V {
        for i in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            vs[0][j][i] = init((i, j), (x, y));
            (vs[0][j][i], trs[0][j][i]) = constraints(t, vs[0][j][i]);
            max_e = trs[0][j][i][0].max(max_e);
        }
    }

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
        boundary: &ghost,
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
        opt: (etaovers, zetaovers, entropy, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: Some(freezeout_energy),
    };

    let e = 2e-1;
    // let e = 2e-1;
    // let e = e.powf(1.0 / r.order as f64);
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_BOTH_2D]; V]; V]; 1],
                   _trs: &[[[[f64; C_MILNE_BOTH_2D]; V]; V]; 1]| {
        let m = vs[0]
            .iter()
            .flat_map(|v| v.iter().map(|v| v[0]))
            .sum::<f64>()
            / (V * V) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order + 1)
    };

    let observables: [Observable<F_BOTH_2D, C_MILNE_BOTH_2D, V, V, 1>; 1] =
        [("momentum_anisotropy", &momentum_anisotropy::<V, V>)];

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
        &err_thr,
        save_raw,
    )
}
