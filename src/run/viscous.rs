use crate::{
    hydro::{
        eos::{hotqcd, wb, EOSs},
        utils::{converge, prepare_trento},
        viscous::viscous2d::{init_from_entropy_density_2d, viscous2d},
        Eos, HydroOutput, C_BOTH_2D, F_BOTH_2D,
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
    init_s: ([[f64; V]; V], usize),
) -> HydroOutput<V, V, F_BOTH_2D, C_BOTH_2D> {
    let (s, i) = init_s;
    let name = format!("InitTrento{}", i);
    let eos = EOSs::WB;
    let (p, dpde, temp): (Eos, Eos, Eos) = match eos {
        EOSs::WB => (&wb::p, &wb::dpde, &wb::T),
        EOSs::HotQCD => (&hotqcd::p, &hotqcd::dpde, &hotqcd::T),
    };
    println!("{}", name);
    let init = init_from_entropy_density_2d(t0, s, p, dpde, temp);
    viscous2d::<V, S>(
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
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let dx = l / V as f64;
    println!("{}", r(0.0).name);
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
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
            // |dt| gauss_legendre_1(Some((0, dt * 0.1))),
            nb_trento,
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
) {
    let trentos = prepare_trento::<V>(nb_trento);
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
        );
    };
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
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
