use crate::{
    boxarray,
    hydro::{utils::eigenvaluesk, Viscosity, C_BOTH_3D, F_BOTH_3D, HBARC},
    solver::{
        context::{BArr, Boundary, Context, Integration, DIM},
        run,
        space::{kt::kt, Dir, Eigenvalues},
        time::{newton::newton, schemes::Scheme},
        utils::{ghost, zeros},
        Constraint, Observable,
    },
};

use crate::hydro::{Eos, Init3D, VOID};

use nalgebra::matrix;

pub fn init_from_entropy_density_3d<'a, const VX: usize, const VY: usize, const VZ: usize>(
    t0: f64,
    s: &'a [[[f64; VX]; VY]; VZ],
    p: Eos<'a>,
    dpde: Eos<'a>,
) -> Box<dyn Fn((usize, usize, usize), (f64, f64, f64)) -> [f64; F_BOTH_3D] + 'a> {
    Box::new(move |(i, j, k), _| {
        let s = s[k][j][i].max(VOID);
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
        move |t, cur @ [tt00, tt01, tt02, tt03, utpi11, utpi12, utpi13, utpi22, utpi23, utpi33, utppi]| {
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
                let ut = (1.0 + ux * ux + uy * uy + uz * uz).sqrt();

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
                let pi33 = utpi33 / ut;
                let pi01 = (ux * pi11 + uy * pi12 + uz*pi13) / ut;
                let pi02 = (ux * pi12 + uy * pi22 + uz*pi23) / ut;
                let pi03 = (ux * pi13 + uy * pi23 + uz*pi33) / ut;
                let pi00 = (ux * pi01 + uy * pi02 + uz*pi03) / ut;
                let pi33 = pi00 - pi01 - pi02 - pi03;

                let m = matrix![
                    pi11, pi12, pi13;
                    pi12, pi22, pi23;
                    pi13, pi23, pi33;
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
                    pi33 * ut,
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
    t: f64,
    [_, _, dpde, ut, _, _, uz, _, _, _, _, _, _, _, _, _, _, _]: [f64; C_BOTH_3D],
) -> f64 {
    eigenvaluesk(dpde, ut, uz) / t
}

pub fn fitutpi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, _pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi]: [f64; C_BOTH_3D],
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
        ut * pi33,
        ut * ppi,
    ]
}
fn fxuxpi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, pi01, _pi02, _pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * ux * ut + pi01),
        t * ((e + pe) * ux * ux + pe + pi11),
        t * ((e + pe) * ux * uy + pi12),
        t * ((e + pe) * ux * uz + pi13),
        ut * pi11,
        ut * pi12,
        ut * pi13,
        ut * pi22,
        ut * pi23,
        ut * pi33,
        ux * ppi,
    ]
}
fn fyuypi(
    t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, pi02, _pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        t * ((e + pe) * uy * ut + pi02),
        t * ((e + pe) * uy * ux + pi12),
        t * ((e + pe) * uy * uy + pe + pi22),
        t * ((e + pe) * uy * uz + pi23),
        ut * pi11,
        ut * pi12,
        ut * pi13,
        ut * pi22,
        ut * pi23,
        ut * pi33,
        uy * ppi,
    ]
}

fn fzuzpi(
    _t: f64,
    [e, pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi]: [f64; C_BOTH_3D],
) -> [f64; F_BOTH_3D] {
    let pe = pe + ppi;
    [
        (e + pe) * uz * ut + pi03,
        (e + pe) * uz * ux + pi13,
        (e + pe) * uz * uy + pi23,
        (e + pe) * uz * uz + pe + pi33,
        ut * pi11,
        ut * pi12,
        ut * pi13,
        ut * pi22,
        ut * pi23,
        ut * pi33,
        uy * ppi,
    ]
}

fn u(
    _t: f64,
    [_e, _pe, _, ut, ux, uy, uz, _pi00, _pi01, _pi02, _pi03, _pi11, _pi12, _pi13, _pi22, _pi23, _pi33, _ppi]: [f64;
        C_BOTH_3D],
) -> [f64; 4] {
    [ut, ux, uy, uz]
}

fn flux<const XY: usize, const Z: usize>(
    [_ov, vs]: [&[[[[f64; F_BOTH_3D]; XY]; XY]; Z]; 2],
    [otrs, trs]: [&[[[[f64; C_BOTH_3D]; XY]; XY]; Z]; 2],
    constraints: Constraint<F_BOTH_3D, C_BOTH_3D>,
    bound: Boundary<F_BOTH_3D, XY, XY, Z>,
    pos: [i32; DIM],
    dx: f64,
    [ot, t]: [f64; 2],
    [_dt, cdt]: [f64; 2],
    &(
        (etas_min, etas_slope, etas_crv),
        (_zetas_max, _zetas_width, _zetas_peak),
        temperature,
        tempcut,
    ): &((f64, f64, f64), (f64, f64, f64), Eos, f64),
) -> [f64; F_BOTH_3D] {
    let tt = t * t;
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
    let (dzf, dzu) = diff(
        (vs, trs),
        bound,
        pos,
        Dir::Z,
        t,
        &fzuzpi,
        &u,
        constraints,
        Eigenvalues::Analytical(&eigenvaluesz),
        pre,
        post,
        dx,
        theta,
    );

    let z = pos[2] as usize;
    let y = pos[1] as usize;
    let x = pos[0] as usize;
    let [_e, _pe, _dpde, out, oux, ouy, ouz, opi00, opi01, opi02, opi03, _opi11, _opi12, _opi13, _opi22, _opi23, _opi33, oppi] =
        otrs[z][y][x];
    let [e, pe, _dpde, ut, ux, uy, uz, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, ppi] =
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
        [0.0, 0.0, 0.0, -tt],
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
    let s = (e + pe) / temp;
    let tc = 0.154; // GeV

    let vgev = gev.max(tc); // viscous temperature [GeV] blocked at tc
    let etaovers = etas_min + etas_slope * (vgev - tc) * (vgev / tc).powf(etas_crv);
    let mut eta = etaovers * s;
    let taupi = 5.0 * eta / (e + pe) + 1e-100; // the 1e-100 is in case etaovers=0

    // let zetaovers = (zetas_max) / (1.0 + ((vmev - zetas_peak) / zetas_width).powi(2));
    // let mut zeta = zetaovers * s;
    // let tauppi = taupi; // use shear relaxation time for bulk

    // bulk pressure is currently turned off
    let mut zeta = 0.0;
    let tauppi = 1.0;

    if gev < tempcut {
        eta = 0.0;
        zeta = 0.0;
    }

    let mut spi = [0.0f64; 11];
    {
        let mut i = 0;
        // pi^munu
        for a in 0..4 {
            for b in a..4 {
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

    let s00 = (e + pe + ppi) * uz * uz + (pe + ppi) + pi33;
    let s03 = (e + pe + ppi) * ut * uz + pi03;
    [
        -dxf[0] - dyf[0] - dzf[0] - dtpi[0] - s00,
        -dxf[1] - dyf[1] - dzf[1] - dtpi[1],
        -dxf[2] - dyf[2] - dzf[2] - dtpi[2],
        -dxf[3] - dyf[3] - dzf[3] - dtpi[3] - s03,
        -dxf[4] - dyf[4] - dzf[4] + spi[4],
        -dxf[5] - dyf[5] - dzf[5] + spi[5],
        -dxf[6] - dyf[6] - dzf[6] + spi[6],
        -dxf[7] - dyf[7] - dzf[7] + spi[7],
        -dxf[8] - dyf[8] - dzf[8] + spi[8],
        -dxf[9] - dyf[9] - dzf[9] + spi[9],
        -dxf[10] - dyf[10] - dzf[10] + spi[10],
    ]
}

// viscous hydro is in Milne coordinates
pub fn viscous3d<const XY: usize, const Z: usize, const S: usize>(
    name: &str,
    maxdt: f64,
    t: f64,
    tend: f64,
    dx: f64,
    r: Scheme<S>,
    p: Eos,
    dpde: Eos,
    temperature: Eos,
    init: Init3D<F_BOTH_3D>,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    shear_temp_cut: f64,
    freezeout_temp: f64,
) -> Option<(
    (BArr<F_BOTH_3D, XY, XY, Z>, BArr<C_BOTH_3D, XY, XY, Z>),
    f64,
    usize,
    usize,
)> {
    let implicit = match r.integration {
        Integration::FixPoint => true,
        Integration::Explicit => false,
    };
    let constraints = gen_constraints(&p, &dpde, temperature, implicit);

    let mut vs: Box<[[[[f64; F_BOTH_3D]; XY]; XY]; Z]> = boxarray(0.0f64);
    let mut trs: Box<[[[[f64; C_BOTH_3D]; XY]; XY]; Z]> = boxarray(0.0f64);
    let names = (
        [
            "tt00", "tt01", "tt02", "tt03", "utpi11", "utpi12", "utpi13", "utpi22", "utpi23",
            "utpi33", "utPi",
        ],
        [
            "e", "pe", "dpde", "ut", "ux", "uy", "uz", "pi00", "pi01", "pi02", "pi03", "pi11",
            "pi12", "pi13", "pi22", "pi23", "pi33", "Pi",
        ],
    );
    let k: Box<[[[[[f64; F_BOTH_3D]; XY]; XY]; Z]; S]> = boxarray(0.0f64);
    let v2 = ((XY - 1) as f64) / 2.0;
    let v2z = ((Z - 1) as f64) / 2.0;
    let mut max_e = 0.0;
    for k in 0..Z {
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
        boundary: &ghost, // TODO use better boundary
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
        opt: (etaovers, zetaovers, temperature, shear_temp_cut),
        p,
        dpde,
        freezeout_energy: Some(freezeout_energy),
    };

    // let e = 2e-3;
    let e = 1e-1;
    let err_thr = |_t: f64,
                   vs: &[[[[f64; F_BOTH_3D]; XY]; XY]; Z],
                   _trs: &[[[[f64; C_BOTH_3D]; XY]; XY]; Z]| {
        let m = vs
            .iter()
            .flat_map(|v| v.iter().flat_map(|v| v.iter().map(|v| v[0])))
            .sum::<f64>()
            / (XY * XY * Z) as f64;
        let k = m / maxdt;
        e * k * (maxdt / dx).powi(r.order)
    };

    let observables: [Observable<F_BOTH_3D, C_BOTH_3D, XY, XY, Z>; 0] = [];

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
    )
}
