use crate::{
    hydro::{
        eos::{conformal_massless, wb, EOSs},
        utils::{converge, prepare_trento_3d},
        viscous::viscous3d::{init_from_energy_density_3d, viscous3d},
        Eos, HydroOutput, C_BOTH_3D, F_BOTH_3D,
    },
    solver::{context::DIM, time::schemes::*, Solver},
    FLOAT,
};

fn hydro3d<const XY: usize, const Z: usize, const S: usize>(
    t0: FLOAT,
    tend: FLOAT,
    dxs: [FLOAT; DIM],
    maxdt: FLOAT,
    etaovers: (FLOAT, FLOAT, FLOAT),
    zetaovers: (FLOAT, FLOAT, FLOAT),
    tempcut: FLOAT,
    freezeout_temp_gev: FLOAT,
    r: Scheme<S>,
    init_s: (&[[[FLOAT; XY]; XY]; Z], usize),
    save_raw: Option<FLOAT>,
) -> HydroOutput<XY, XY, Z, F_BOTH_3D, C_BOTH_3D> {
    let (s, i) = init_s;
    let name = ("InitTrento", i);
    let eos = EOSs::WB;
    // let eos = EOSs::HotQCD;
    // let eos = EOSs::HotQCDLog;
    let (p, dpde, entropy, temp): (Eos, Eos, Eos, Eos) = match eos {
        EOSs::ConformalMassless => (
            &conformal_massless::p,
            &conformal_massless::dpde,
            &conformal_massless::s,
            &conformal_massless::T,
        ),
        EOSs::WB => (&wb::p, &wb::dpde, &wb::s, &wb::T),
    };
    println!("{}{}", name.0, name.1);
    let init = init_from_energy_density_3d(t0, s, p, dpde, entropy, temp);
    viscous3d::<XY, Z, S>(
        &name,
        maxdt,
        t0,
        tend,
        dxs,
        r,
        p,
        dpde,
        entropy,
        temp,
        &init,
        etaovers,
        zetaovers,
        tempcut,
        freezeout_temp_gev,
        save_raw,
    )
}

pub fn run_convergence_3d<const XY: usize, const Z: usize, const S: usize>(
    t0: FLOAT,
    tend: FLOAT,
    l: FLOAT,
    etas_len: FLOAT,
    dtmin: FLOAT,
    dtmax: FLOAT,
    etaovers: (FLOAT, FLOAT, FLOAT),
    zetaovers: (FLOAT, FLOAT, FLOAT),
    tempcut: FLOAT,
    freezeout_temp_mev: FLOAT,
    r: impl Fn(FLOAT) -> Scheme<S>,
    (nb_trento, first_trento): (usize, usize),
    save_raw: Option<FLOAT>,
) {
    let (trentos, dxs) = prepare_trento_3d::<XY, Z>(nb_trento, first_trento);
    let dx = l / XY as FLOAT;
    let detas = etas_len / Z as FLOAT;
    let dxs = dxs.unwrap_or([dx, dx, detas]);
    println!("{}", r(0.0).name);
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), first_trento + i);
        converge(dtmax, dtmin, |dt| {
            hydro3d::<XY, Z, S>(
                t0,
                tend,
                dxs,
                dt,
                etaovers,
                zetaovers,
                tempcut,
                freezeout_temp_mev,
                r(dt),
                trento,
                save_raw,
            )
        });
    }
}
pub fn run_3d<const XY: usize, const Z: usize>(
    solver: Solver,
    t0: FLOAT,
    tend: FLOAT,
    l: FLOAT,
    etas_len: FLOAT,
    dtmin: FLOAT,
    dtmax: FLOAT,
    etaovers: (FLOAT, FLOAT, FLOAT),
    zetaovers: (FLOAT, FLOAT, FLOAT),
    tempcut: FLOAT,
    freezeout_temp_gev: FLOAT,
    nf_trento: (usize, usize),
    save_raw: Option<FLOAT>,
) {
    let do_gl1 = || {
        run_convergence_3d::<XY, Z, 1>(
            t0,
            tend,
            l,
            etas_len,
            dtmin,
            dtmax,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            |_| gauss_legendre_1(),
            nf_trento,
            save_raw,
        )
    };
    let do_heun = || {
        run_convergence_3d::<XY, Z, 2>(
            t0,
            tend,
            l,
            etas_len,
            dtmin,
            dtmax,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            |_| heun(),
            nf_trento,
            save_raw,
        )
    };
    match solver {
        Solver::Both => {
            do_gl1();
            do_heun();
        }
        Solver::Implicit => {
            do_gl1();
        }
        Solver::Explicit => {
            do_heun();
        }
    }
}
pub fn run_trento_3d<const XY: usize, const Z: usize>(
    solver: Solver,
    t0: FLOAT,
    tend: FLOAT,
    l: FLOAT,
    etas_len: FLOAT,
    dt: FLOAT,
    etaovers: (FLOAT, FLOAT, FLOAT),
    zetaovers: (FLOAT, FLOAT, FLOAT),
    tempcut: FLOAT,
    freezeout_temp_gev: FLOAT,
    (nb_trento, first_trento): (usize, usize),
    save_raw: Option<FLOAT>,
) {
    let (trentos, dxs) = prepare_trento_3d::<XY, Z>(nb_trento, first_trento);
    let gl1 = gauss_legendre_1();
    let heun = heun();
    let dx = l / XY as FLOAT;
    let detas = etas_len / Z as FLOAT;
    let dxs = dxs.unwrap_or([dx, dx, detas]);
    let do_gl1 = |trento| {
        hydro3d::<XY, Z, 1>(
            t0,
            tend,
            dxs,
            dt,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            gl1,
            trento,
            save_raw,
        );
    };
    let do_heun = |trento| {
        hydro3d::<XY, Z, 2>(
            t0,
            tend,
            dxs,
            dt,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            heun,
            trento,
            save_raw,
        );
    };
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), first_trento + i);
        match solver {
            Solver::Both => {
                do_gl1(trento);
                do_heun(trento);
            }
            Solver::Implicit => {
                do_gl1(trento);
            }
            Solver::Explicit => {
                do_heun(trento);
            }
        }
    }
}
