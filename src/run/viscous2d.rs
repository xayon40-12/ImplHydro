use crate::{
    hydro::{
        eos::{conformal_massless, hotqcd, wb, EOSs},
        utils::{converge, prepare_trento_2d},
        viscous::viscous2d::{init_from_entropy_density_2d, viscous2d},
        Eos, HydroOutput, C_MILNE_BOTH_2D, F_BOTH_2D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_gev: f64,
    r: Scheme<S>,
    init_s: (&[[f64; V]; V], usize),
    save_raw: bool,
) -> HydroOutput<V, V, 1, F_BOTH_2D, C_MILNE_BOTH_2D> {
    let (s, i) = init_s;
    let name = ("InitTrento", i);
    // let eos = EOSs::ConformalMassless;
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
        EOSs::HotQCD => (&hotqcd::p, &hotqcd::dpde, &hotqcd::s, &hotqcd::T),
        EOSs::HotQCDLog => (
            &hotqcd::log::p,
            &hotqcd::log::dpde,
            &hotqcd::log::s,
            &hotqcd::log::T,
        ),
    };
    println!("{}{}", name.0, name.1);
    let init = init_from_entropy_density_2d(t0, s, p, dpde, entropy);
    viscous2d::<V, S>(
        &name,
        maxdt,
        t0,
        tend,
        dx,
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

pub fn run_convergence_2d<const V: usize, const S: usize>(
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
    save_raw: bool,
) {
    let trentos = prepare_trento_2d::<V>(nb_trento);
    let dx = l / V as f64;
    println!("{}", r(0.0).name);
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), i);
        converge(dtmax, dtmin, |dt| {
            hydro2d::<V, S>(
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
                save_raw,
            )
        });
    }
}
pub fn run_2d<const V: usize>(
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
    save_raw: bool,
) {
    let do_gl1 = || {
        run_convergence_2d::<V, 1>(
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
            save_raw,
        )
    };
    let do_heun = || {
        run_convergence_2d::<V, 2>(
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
pub fn run_trento_2d<const V: usize>(
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
    save_raw: bool,
) {
    let trentos = prepare_trento_2d::<V>(nb_trento);
    let gl1 = gauss_legendre_1();
    let heun = heun();
    let dx = l / V as f64;
    let do_gl1 = |trento| {
        hydro2d::<V, 1>(
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
            save_raw,
        );
    };
    let do_heun = |trento| {
        hydro2d::<V, 2>(
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
            save_raw,
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
