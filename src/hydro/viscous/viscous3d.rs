use crate::{
    hydro::{utils::eigenvaluesk, Viscosity, C_BOTH_3D, F_BOTH_3D, HBARC},
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

use crate::hydro::{Eos, Init3D, VOID};

use nalgebra::matrix;

pub fn init_from_entropy_density_3d<'a, const VX: usize, const VY: usize, const VZ: usize>(
    t0: f64,
    s: &'a [[[f64; VX]; VY]; VZ],
    p: Eos<'a>,
    dpde: Eos<'a>,
    entropy: Eos<'a>,
) -> Box<dyn Fn((usize, usize, usize), (f64, f64, f64)) -> [f64; F_BOTH_3D] + 'a> {
    Box::new(move |(i, j, k), _| {
        let s = s[k][j][i];
        let e = newton(1e-10, s, |e| entropy(e) - s, |e| e.max(0.0).min(1e10)).max(VOID);
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
        ];
        fitutpi(t0, vars)
    })
}

fn gen_constraints<'a>(
    p: Eos<'a>,
    dpde: Eos<'a>,
    _temp: Eos<'a>,
    _implicit: bool,
) -> Box<dyn Fn(f64, [f64; F_BOTH_3D]) -> ([f64; F_BOTH_3D], [f64; C_BOTH_3D]) + 'a + Sync> {
    Box::new(
        move |t, cur @ [tt00, tt01, tt02, tt03, utpi11, utpi12, utpi13, utpi22, utpi23, utppi]| {
            let t00 = tt00 / t;
            let t01 = tt01 / t;
            let t02 = tt02 / t;
            let t03 = tt03 / t;
            let m = (t01 * t01 + t02 * t02 + t03 * t03).sqrt();
            let t00 = t00.max(m * (1.0 + 1e-15));

            let sv = |v: f64| m / (t00 + p((t00 - m * v).max(VOID)));
            let cv = |v: f64| v.max(0.0).min(1.0);
            let gamma_max2 = 20.0f64.powi(2);
            let v_max = (1.0 - 1.0 / gamma_max2).sqrt();
            let v = newton(1e-10, 0.5, |v| sv(v) - v, cv).min(v_max);

            let e = (t00 - m * v).max(VOID).min(1e10);

            {
                let g = (1.0 - v * v).sqrt();

                let pe = p(e);
                let ux = g * t01 / (e + pe);
                let uy = g * t02 / (e + pe);
                let uz = g * t03 / (e + pe);
                let ux2 = ux * ux;
                let uy2 = uy * uy;
                let uz2 = uz * uz;
                let ut2 = 1.0 + ux2 + uy2 + uz2;
                let ut = ut2.sqrt();

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
                let pi13 = utpi13 / ut;
                let pi22 = utpi22 / ut;
                let pi23 = utpi23 / ut;
                let pi33 = (2.0
                    * (ux * uy * pi12
                        + ux * uz * pi13
                        + uy * ux * pi12
                        + uy * uz * pi23
                        + uz * ux * pi13
                        + uz * uy * pi23)
                    + (ux2 - ut2) * pi11
                    + (uy2 - ut2) * pi22)
                    / (ut2 - uz2);
                let pi01 = (ux * pi11 + uy * pi12 + uz * pi13) / ut;
                let pi02 = (ux * pi12 + uy * pi22 + uz * pi23) / ut;
                let pi03 = (ux * pi13 + uy * pi23 + uz * pi33) / ut;
                let pi00 = (ux * pi01 + uy * pi02 + uz * pi03) / ut;
                let pi33 = pi00 - pi11 - pi22;

                let m = matrix![
                    pi00, pi01, pi02, pi03;
                    pi01, pi11, pi12, pi13;
                    pi02, pi12, pi22, pi23;
                    pi03, pi13, pi23, pi33;
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
                let pi13 = r * pi13;
                let pi22 = r * pi22;
                let pi23 = r * pi23;
                let pi33 = r * pi33;
                let pi01 = r * pi01;
                let pi02 = r * pi02;
                let pi03 = r * pi03;
                let pi00 = r * pi00;

                let vs = [
                    t00 * t,
                    t01 * t,
                    t02 * t,
                    t03 * t,
                    pi11 * ut,
                    pi12 * ut,
                    pi13 * ut,
                    pi22 * ut,
                    pi23 * ut,
                    ppi * ut,
                ];

                let trans = [
                    e,
                    pe,
                    dpde(e),
                    ut,
                    ux,
                    uy,
                    uz,
                    pi00,
                    pi01,
                    pi02,
                    pi03,
                    pi11,
                    pi12,
                    pi13,
                    pi22,
                    pi23,
                    pi33,
                    ppi,
                ];

                if vs.iter().any(|v| v.is_nan()) || trans.iter().any(|v| v.is_nan()) {
                    panic!(
                        "\n\nNaN in constraint\n{:?}\n{} {}\n{:?}\n{:?}\n\n",
                        cur, g, v, vs, trans
                    );
                }

                (vs, trans)
            }
        },
    )
}

fn eigenvaluesx(
    _t: f64,
    [_, _, dpde, ut, ux, _, _, _, _, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_3D],
) -> f64 {
    eigenvaluesk(dpde, ut, ux)
}
fn eigenvaluesy(
    _t: f64,
    [_, _, dpde, ut, _, uy, _, _, _, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_3D],
) -> f64 {
    eigenvaluesk(dpde, ut, uy)
}
fn eigenvaluesz(
    _t: f64,
    [_, _, dpde, ut, _, _, uz, _, _, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_3D],
) -> f64 {
    eigenvaluesk(dpde, ut, uz)
}

pub fn fitutpi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, _pi03, pi11, pi12, pi13, pi22, pi23, _pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    [
        t * ((e + pe) * ut * ut - pe), // no pi in Ttt because it is substracted in the flux
        t * ((e + pe) * ut * ux),
        t * ((e + pe) * ut * uy),
        t * ((e + pe) * ut * uz),
        ut * pi11,
        ut * pi12,
        ut * pi13,
        ut * pi22,
        ut * pi23,
        ut * ppi,
    ]
}
fn fxuxpi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, pi01, _pi02, _pi03, pi11, pi12, pi13, pi22, pi23, _pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * ux * ut + pi01),
        t * ((e + pe) * ux * ux + pe + pi11),
        t * ((e + pe) * ux * uy + pi12),
        t * ((e + pe) * ux * uz + pi13),
        ux * pi11,
        ux * pi12,
        ux * pi13,
        ux * pi22,
        ux * pi23,
        ux * ppi,
    ]
}
fn fyuypi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, pi02, _pi03, pi11, pi12, pi13, pi22, pi23, _pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * uy * ut + pi02),
        t * ((e + pe) * uy * ux + pi12),
        t * ((e + pe) * uy * uy + pe + pi22),
        t * ((e + pe) * uy * uz + pi23),
        uy * pi11,
        uy * pi12,
        uy * pi13,
        uy * pi22,
        uy * pi23,
        uy * ppi,
    ]
}

fn fzuzpi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * uz * ut + pi03),
        t * ((e + pe) * uz * ux + pi13),
        t * ((e + pe) * uz * uy + pi23),
        t * ((e + pe) * uz * uz + pe + pi33),
        uz * pi11,
        uz * pi12,
        uz * pi13,
        uz * pi22,
        uz * pi23,
        uz * ppi,
    ]
}

fn u(
    _t: f64,
    [_e, _pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, _pi03, _pi11, _pi12, _pi13, _pi22, _pi23, _pi33, _ppi]: [f64;
        C_BOTH_3D],
) -> [f64; 4] {
    [ut, ux, uy, uz]
}

fn flux<const XY: usize, const VZ: usize>(
    k: &Arr<F_BOTH_3D, XY, XY, VZ>,
    [_ov, vs]: [&Arr<F_BOTH_3D, XY, XY, VZ>; 2],
    [otrs, trs]: [&Arr<C_BOTH_3D, XY, XY, VZ>; 2],
    constraints: Constraint<F_BOTH_3D, C_BOTH_3D>,
    bound: Boundary<F_BOTH_3D, XY, XY, VZ>,
    pos: [i32; DIM],
    dx: f64,
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    &(
        (etas_min, etas_slope, etas_crv),
        (zetas_max, zetas_width, zetas_peak),
        entropy,
        temperature,
        tempcut,
    ): &((f64, f64, f64), (f64, f64, f64), Eos, Eos, f64),
) -> [f64; F_BOTH_3D] {
    let theta = 1.1;

    let pre = &|_t: f64, mut vs: [f64; F_BOTH_3D]| {
        let t00 = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let t03 = vs[3];
        let k = t01 * t01 + t02 * t02 + t03 * t03;
        let m = (t00 * t00 - k).sqrt();
        vs[0] = m;
        vs
    };
    let post = &|_t: f64, mut vs: [f64; F_BOTH_3D]| {
        let m = vs[0];
        let t01 = vs[1];
        let t02 = vs[2];
        let t03 = vs[3];
        let k = t01 * t01 + t02 * t02 + t03 * t03;
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
        Z(FluxInfo {
            flux: &fzuzpi,
            secondary: &u,
            eigenvalues: Eigenvalues::Analytical(&eigenvaluesz),
        }),
    ];
    let [(dxf, dxu), (dyf, dyu), (mut dzf, mut dzu)] = diff(
        (vs, trs),
        bound,
        pos,
        t,
        flux_infos,
        constraints,
        pre,
        post,
        dx,
        theta,
    );
    // in this formalism `dz` becomes `(1/t)*dz`
    for i in 0..dzf.len() {
        dzf[i] /= t;
    }
    for i in 0..dzu.len() {
        dzu[i] /= t;
    }

    let z = pos[2] as usize;
    let y = pos[1] as usize;
    let x = pos[0] as usize;
    let [ktt00, ktt01, ktt02, ktt03, kutpi11, kutpi12, kutpi13, kutpi22, kutpi23, kutppi] =
        k[z][y][x];
    let [tt00, tt01, tt02, tt03, _utpi11, _utpi12, _utpi13, _utpi22, _utpi23, _utppi] = vs[z][y][x];
    let [_oe, _ope, _odpde, out, oux, ouy, ouz, opi00, opi01, opi02, opi03, _opi11, _opi12, _opi13, _opi22, _opi23, _opi33, oppi] =
        otrs[z][y][x];
    let [e, pe, dpde, ut, ux, uy, uz, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi] =
        trs[z][y][x];
    let pimn = [
        [pi00, pi01, pi02, pi03],
        [pi01, pi11, pi12, pi13],
        [pi02, pi12, pi22, pi23],
        [pi03, pi13, pi23, pi33],
    ];

    let g = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
    ];
    let u = [ut, ux, uy, uz];
    let mut delta = [[0.0f64; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            delta[a][b] = g[a][b] - u[a] * u[b];
        }
    }
    let ou = [out, oux, ouy, ouz];
    let mut odelta = [[0.0f64; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            odelta[a][b] = g[a][b] - ou[a] * ou[b];
        }
    }

    let (dtu, dtpi) = {
        if EXACT {
            let ut2 = ut * ut;
            let ux2 = ux * ux;
            let uy2 = uy * uy;
            let uz2 = uz * uz;
            let d0 = dpde * ut2 - dpde + ut2;
            let dp0 = dpde * ut2 + dpde + ut2;
            let d12 = d0 * ut2 - dp0 * ux2 - dp0 * uy2;
            let d13 = d0 * ut2 - dp0 * ux2 - dp0 * uz2;
            let d23 = d0 * ut2 - dp0 * uy2 - dp0 * uz2;
            let d123 = d0 * ut2 - dp0 * ux2 - dp0 * uy2 - dp0 * uz2;
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
            //     ⎡  2     2     2     2                                                            ⎤
            //     ⎢u₀  + u₁  + u₂  + u₃       -2⋅u₀⋅u₁            -2⋅u₀⋅u₂            -2⋅u₀⋅u₃      ⎥
            //     ⎢─────────────────────      ─────────           ─────────           ─────────     ⎥
            //     ⎢         D₁₂₃                 D₁₂₃                D₁₂₃                D₁₂₃       ⎥
            //     ⎢                                                                                 ⎥
            //     ⎢    2                                                                            ⎥
            //     ⎢ -u₀ ⋅u₁⋅(dpde + 1)           D₂₃              DP₀⋅u₁⋅u₂           DP₀⋅u₁⋅u₃     ⎥
            //     ⎢ ───────────────────   ──────────────────  ──────────────────  ──────────────────⎥
            //     ⎢   D₁₂₃⋅(e + p(e))     D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))⎥
            // m = ⎢                                                                                 ⎥
            //     ⎢    2                                                                            ⎥
            //     ⎢ -u₀ ⋅u₂⋅(dpde + 1)        DP₀⋅u₁⋅u₂              D₁₃              DP₀⋅u₂⋅u₃     ⎥
            //     ⎢ ───────────────────   ──────────────────  ──────────────────  ──────────────────⎥
            //     ⎢   D₁₂₃⋅(e + p(e))     D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))⎥
            //     ⎢                                                                                 ⎥
            //     ⎢    2                                                                            ⎥
            //     ⎢ -u₀ ⋅u₃⋅(dpde + 1)        DP₀⋅u₁⋅u₃           DP₀⋅u₂⋅u₃              D₁₂        ⎥
            //     ⎢ ───────────────────   ──────────────────  ──────────────────  ──────────────────⎥
            //     ⎣   D₁₂₃⋅(e + p(e))     D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))  D₁₂₃⋅u₀⋅(e + p(e))⎦

            let d123m = [
                [
                    ut2 + ux2 + uy2 + uz2,
                    -2.0 * ut * ux,
                    -2.0 * ut * uy,
                    -2.0 * ut * uz,
                ],
                [
                    -ut2 * ux * dpde1 / ep,
                    d23 / u0ep,
                    dp0 * ux * uy / u0ep,
                    dp0 * ux * uz / u0ep,
                ],
                [
                    -ut2 * uy * dpde1 / ep,
                    dp0 * ux * uy / u0ep,
                    d13 / u0ep,
                    dp0 * uy * uz / u0ep,
                ],
                [
                    -ut2 * uz * dpde1 / ep,
                    dp0 * ux * uz / u0ep,
                    dp0 * uy * uz / u0ep,
                    d12 / u0ep,
                ],
            ];
            // d_t (tT^{t\nu}) = T^{t\nu} + t d_t T^{t\nu}
            let k = [
                (ktt00 - tt00 / t) / t,
                (ktt01 - tt01 / t) / t,
                (ktt02 - tt02 / t) / t,
                (ktt03 - tt03 / t) / t,
            ];
            let mut v = [0.0f64; 4];
            for j in 0..4 {
                for i in 0..4 {
                    v[j] += d123m[j][i] / d123 * k[i];
                }
            }
            let _dedt = v[0];
            let dtu = [(ux * v[1] + uy * v[2] + uz * v[3]) / ut, v[1], v[2], v[3]];
            let udtu = [ut * dtu[0], ux * dtu[1], uy * dtu[2], uz * dtu[3]];
            // d_t(ut*Pi) = Pi*d_t ut + ut*d_t Pi
            let dtppi = (kutppi - ppi * dtu[0]) / ut;
            // d_t(ut*pi^{t\nu}) = pi^{t\nu}*d_t ut + ut*d_t pi^{t\nu}
            let dtpi11 = (kutpi11 - pimn[1][1] * dtu[0]) / ut;
            let dtpi12 = (kutpi12 - pimn[1][2] * dtu[0]) / ut;
            let dtpi13 = (kutpi13 - pimn[1][3] * dtu[0]) / ut;
            let dtpi22 = (kutpi22 - pimn[2][2] * dtu[0]) / ut;
            let dtpi23 = (kutpi23 - pimn[2][3] * dtu[0]) / ut;
            // let dtpi33 = (kutpi33 - pimn[3][3] * dtu[0]) / ut;
            let dtcross = 2.0
                * ((ux * uy * dtpi12
                    + ux * uz * dtpi13
                    + uy * ux * dtpi12
                    + uy * uz * dtpi23
                    + uz * ux * dtpi13
                    + uz * uy * dtpi23)
                    + (dtu[1] * uy * pi12
                        + dtu[1] * uz * pi13
                        + dtu[2] * ux * pi12
                        + dtu[2] * uz * pi23
                        + dtu[3] * ux * pi13
                        + dtu[3] * uy * pi23)
                    + (ux * dtu[2] * pi12
                        + ux * dtu[3] * pi13
                        + uy * dtu[1] * pi12
                        + uy * dtu[3] * pi23
                        + uz * dtu[1] * pi13
                        + uz * dtu[2] * pi23)
                    + (udtu[1] - udtu[0]) * pi11
                    + (udtu[2] - udtu[0]) * pi22)
                + (ux2 - ut2) * dtpi11
                + (uy2 - ut2) * dtpi22;
            let dtpi33 = (dtcross - pimn[3][3] * 2.0 * (udtu[0] - udtu[3])) / (ut2 - uz2);

            // utpi01 = uxpi11 + uypi12 + uzpi13
            // kutpi01 = uxdtpi11 + pi11dtux + uydtpi12 + pi12dtuy + uzdtpi13 + pi13dtuz
            let kutpi01 = ux * dtpi11
                + pi11 * dtu[1]
                + uy * dtpi12
                + pi12 * dtu[2]
                + uz * dtpi13
                + pi13 * dtu[3];
            let dtpi01 = (kutpi01 - pimn[1][0] * dtu[0]) / ut;
            // utpi02 = uxpi21 + uypi22 + uzpi23
            let kutpi02 = ux * dtpi12
                + pi12 * dtu[1]
                + uy * dtpi22
                + pi22 * dtu[2]
                + uz * dtpi23
                + pi23 * dtu[3];
            let dtpi02 = (kutpi02 - pimn[2][0] * dtu[0]) / ut;
            // utpi03 = uxpi13 + uypi23 + uzpi33
            let kutpi03 = ux * dtpi13
                + pi13 * dtu[1]
                + uy * dtpi23
                + pi23 * dtu[2]
                + uz * dtpi33
                + pi33 * dtu[3];
            let dtpi03 = (kutpi03 - pimn[3][0] * dtu[0]) / ut;
            // utpi00 = uxpi10 + uypi20
            let kutpi00 = ux * dtpi01
                + pi01 * dtu[1]
                + uy * dtpi02
                + pi02 * dtu[2]
                + uz * dtpi03
                + pi03 * dtu[3];
            let dtpi00 = (kutpi00 - pimn[0][0] * dtu[0]) / ut;
            let dttpi00 = pi00 + t * dtpi00;
            let dttpi01 = pi01 + t * dtpi01;
            let dttpi02 = pi02 + t * dtpi02;
            let dttpi03 = pi03 + t * dtpi03;
            let mut dttppideltat = [0.0f64; 4];
            for i in 0..4 {
                dttppideltat[i] = ppi * delta[0][i]
                    + t * (delta[0][i] * dtppi - ppi * (u[i] * dtu[0] + u[0] * dtu[i]));
            }
            let dtpi = [
                dttpi00 - dttppideltat[0],
                dttpi01 - dttppideltat[1],
                dttpi02 - dttppideltat[2],
                dttpi03 - dttppideltat[3],
            ];

            (dtu, dtpi)
        } else {
            let dtu = [
                (ut - out) / cdt,
                (ux - oux) / cdt,
                (uy - ouy) / cdt,
                (uz - ouz) / cdt,
            ];
            let dtpi = [
                (t * (pi00 - ppi * delta[0][0]) - ot * (opi00 - oppi * odelta[0][0])) / cdt,
                (t * (pi01 - ppi * delta[0][1]) - ot * (opi01 - oppi * odelta[0][1])) / cdt,
                (t * (pi02 - ppi * delta[0][2]) - ot * (opi02 - oppi * odelta[0][2])) / cdt,
                (t * (pi03 - ppi * delta[0][3]) - ot * (opi03 - oppi * odelta[0][3])) / cdt,
            ];
            (dtu, dtpi)
        }
    };
    let du = [dtu, dxu, dyu, dzu];
    let mut ddu = [0.0f64; 4];
    let mut dcuc = u[0] / t;
    for j in 0..4 {
        dcuc += du[j][j];
        for i in 0..4 {
            ddu[j] += u[i] * du[i][j];
        }
    }

    let temp = temperature(e);
    let gev = temp * HBARC;
    let s = entropy(e); // (e + pe) / temp;
    let tc = 0.154; // GeV

    let vgev = gev.max(tc); // viscous temperature [GeV] blocked at tc
    let etaovers = etas_min + etas_slope * (vgev - tc) * (vgev / tc).powf(etas_crv);
    let mut eta = etaovers * s;
    let taupi = 5.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0

    let zetaovers = (zetas_max) / (1.0 + ((vgev - zetas_peak) / zetas_width).powi(2));
    let mut zeta = zetaovers * s;
    let tauppi = taupi; // use shear relaxation time for bulk

    if gev < tempcut {
        eta = 0.0;
        zeta = 0.0;
    }

    // t * I^{\mu\nu}_{\pi,G}
    let t_ipi_g = [
        [
            2.0 * uz * pimn[0][3],
            uz * pimn[3][1],
            uz * pimn[3][2],
            uz * (pimn[0][0] + pimn[3][3]),
        ],
        [uz * pimn[3][1], 0.0, 0.0, uz * pimn[0][1]],
        [uz * pimn[3][2], 0.0, 0.0, uz * pimn[0][2]],
        [
            uz * (pimn[0][0] + pimn[3][3]),
            uz * pimn[0][1],
            uz * pimn[0][2],
            2.0 * uz * pimn[0][3],
        ],
    ];
    let mut spi = [0.0f64; 11];
    {
        let mut i = 0;
        // pi^munu
        for a in 0..4 {
            for b in a..4 {
                if !(a == 3 && b == 3) {
                    let pi_ns = eta
                        * ((0..4)
                            .map(|i| delta[a][i] * du[i][b] + delta[b][i] * du[i][a])
                            .sum::<f64>()
                            - 2.0 / 3.0 * delta[a][b] * dcuc);
                    let ipi = -1.0 / 3.0 * pimn[a][b] * dcuc
                        - pimn[a][b] * u[0] / t
                        - (0..4)
                            .map(|i| (u[a] * pimn[b][i] + u[b] * pimn[a][i]) * ddu[i] * g[i][i])
                            .sum::<f64>()
                        - t_ipi_g[a][b] / t;

                    spi[i] = -(pimn[a][b] - pi_ns) / taupi + ipi;
                    i += 1;
                }
            }
        }
        // Pi
        let ppi_ns = -zeta * dcuc;
        let ippi = -1.0 / 3.0 * ppi * dcuc - ppi * u[0] / t;
        spi[i] = -(ppi - ppi_ns) / tauppi + ippi;
    }

    let smunu = |i: usize, j: usize| e * u[i] * u[j] - (pe + ppi) * delta[i][j] + pimn[i][j];
    [
        -dxf[0] - dyf[0] - dzf[0] - dtpi[0] - smunu(3, 3),
        -dxf[1] - dyf[1] - dzf[1] - dtpi[1],
        -dxf[2] - dyf[2] - dzf[2] - dtpi[2],
        -dxf[3] - dyf[3] - dzf[3] - dtpi[3] - smunu(0, 3),
        -dxf[4] - dyf[4] - dzf[4] + spi[4],
        -dxf[5] - dyf[5] - dzf[5] + spi[5],
        -dxf[6] - dyf[6] - dzf[6] + spi[6],
        -dxf[7] - dyf[7] - dzf[7] + spi[7],
        -dxf[8] - dyf[8] - dzf[8] + spi[8],
        -dxf[9] - dyf[9] - dzf[9] + spi[9],
    ]
}

// viscous hydro is in Milne coordinates
pub fn viscous3d<const XY: usize, const VZ: usize, const S: usize>(
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
    init: Init3D<F_BOTH_3D>,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    shear_temp_cut: f64,
    freezeout_temp: f64,
    save_raw: Option<f64>,
) -> Option<(
    (BArr<F_BOTH_3D, XY, XY, VZ>, BArr<C_BOTH_3D, XY, XY, VZ>),
    f64,
    usize,
    usize,
)> {
    let implicit = match r.integration {
        Integration::FixPoint => true,
        Integration::Explicit => false,
    };
    let constraints = gen_constraints(&p, &dpde, temperature, implicit);

    let mut vs: Box<[[[[f64; F_BOTH_3D]; XY]; XY]; VZ]> = boxarray(0.0);
    let mut trs: Box<[[[[f64; C_BOTH_3D]; XY]; XY]; VZ]> = boxarray(0.0);
    let names = (
        [
            "tt00", "tt01", "tt02", "tt03", "utpi11", "utpi12", "utpi13", "utpi22", "utpi23",
            "utPi",
        ],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "uz", "pi00", "pi01", "pi02", "pi03", "pi11",
            "pi12", "pi13", "pi22", "pi23", "pi33", "Pi",
        ],
    );
    let k: Box<[[[[[f64; F_BOTH_3D]; XY]; XY]; VZ]; S]> = boxarray(0.0);
    let v2 = ((XY - 1) as f64) / 2.0;
    let v2z = ((VZ - 1) as f64) / 2.0;
    let mut max_e = 0.0;
    for k in 0..VZ {
        for j in 0..XY {
            for i in 0..XY {
                let x = (i as f64 - v2) * dx;
                let y = (j as f64 - v2) * dx;
                let z = (k as f64 - v2z) * dx;
                vs[k][j][i] = init((i, j, k), (x, y, z));
                (vs[k][j][i], trs[k][j][i]) = constraints(t, vs[k][j][i]);
                max_e = trs[k][j][i][0].max(max_e);
            }
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
        opt: (etaovers, zetaovers, entropy, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: Some(freezeout_energy),
    };

    let e = 5e-2;
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_BOTH_3D]; XY]; XY]; VZ],
                   _trs: &[[[[f64; C_BOTH_3D]; XY]; XY]; VZ]| {
        let m = vs
            .iter()
            .flat_map(|v| v.iter().flat_map(|v| v.iter().map(|v| v[0])))
            .sum::<f64>()
            / (XY * XY * VZ) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order + 1)
    };

    let observables: [Observable<F_BOTH_3D, C_BOTH_3D, XY, XY, VZ>; 0] = [];

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
