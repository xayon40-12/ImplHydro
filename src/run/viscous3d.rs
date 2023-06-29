use crate::{
    hydro::{
        eos::{hotqcd, wb, EOSs},
        utils::{converge, prepare_trento_3d},
        viscous::viscous3d::{init_from_entropy_density_3d, viscous3d},
        Eos, HydroOutput, C_BOTH_3D, F_BOTH_3D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro3d<const XY: usize, const Z: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_gev: f64,
    r: Scheme<S>,
    init_s: (&[[[f64; XY]; XY]; Z], usize),
) -> HydroOutput<XY, XY, Z, F_BOTH_3D, C_BOTH_3D> {
    let (s, i) = init_s;
    let name = format!("InitTrento{}", i);
    let eos = EOSs::WB;
    // let eos = EOSs::HotQCD;
    let (p, dpde, temp): (Eos, Eos, Eos) = match eos {
        EOSs::WB => (&wb::p, &wb::dpde, &wb::T),
        EOSs::HotQCD => (&hotqcd::p, &hotqcd::dpde, &hotqcd::T),
    };
    println!("{}", name);
    let init = init_from_entropy_density_3d(t0, s, p, dpde);
    viscous3d::<XY, Z, S>(
        &name,
        maxdt,
        t0,
        tend,
        dx,
        r,
        p,
        dpde,
        temp,
        &init,
        etaovers,
        zetaovers,
        tempcut,
        freezeout_temp_gev,
    )
}

pub fn run_convergence_3d<const XY: usize, const Z: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_mev: f64,
    r: impl Fn(f64) -> Scheme<S>,
    nb_trento: usize,
) {
    let trentos = prepare_trento_3d::<XY, Z>(nb_trento);
    let dx = l / XY as f64;
    println!("{}", r(0.0).name);
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), i);
        converge(dtmax, dtmin, |dt| {
            hydro3d::<XY, Z, S>(
                t0,
                tend,
                dx,
                dt,
                etaovers,
                zetaovers,
                tempcut,
                freezeout_temp_mev,
                r(dt),
                trento,
            )
        });
    }
}
pub fn run_3d<const XY: usize, const Z: usize>(
    solver: Solver,
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_gev: f64,
    nb_trento: usize,
) {
    let do_gl1 = || {
        run_convergence_3d::<XY, Z, 1>(
            t0,
            tend,
            l,
            dtmin,
            dtmax,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            |_| gauss_legendre_1(),
            nb_trento,
        )
    };
    let do_heun = || {
        run_convergence_3d::<XY, Z, 2>(
            t0,
            tend,
            l,
            dtmin,
            dtmax,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            |_| heun(),
            nb_trento,
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
    t0: f64,
    tend: f64,
    l: f64,
    dt: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_gev: f64,
    nb_trento: usize,
) {
    let trentos = prepare_trento_3d::<XY, Z>(nb_trento);
    let gl1 = gauss_legendre_1();
    let heun = heun();
    let dx = l / XY as f64;
    let do_gl1 = |trento| {
        hydro3d::<XY, Z, 1>(
            t0,
            tend,
            dx,
            dt,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            gl1,
            trento,
        );
    };
    let do_heun = |trento| {
        hydro3d::<XY, Z, 2>(
            t0,
            tend,
            dx,
            dt,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_gev,
            heun,
            trento,
        );
    };
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), i);
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
