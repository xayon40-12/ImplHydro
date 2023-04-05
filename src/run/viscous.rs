use crate::{
    hydro::{
        eos::{hotqcd, wb, EOSs},
        utils::{converge, prepare_trento_freestream},
        viscous::init_from_freestream_2d,
        viscous::viscous2d::viscous2d,
        Eos, HydroOutput, C_BOTH_2D, FREESTREAM_2D, F_BOTH_2D,
    },
    solver::time::schemes::*,
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    er: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_mev: f64,
    r: Scheme<S>,
    init_s: ([[[f64; FREESTREAM_2D]; V]; V], usize),
) -> HydroOutput<V, V, F_BOTH_2D, C_BOTH_2D> {
    let (s, i) = init_s;
    let name = format!("InitTrento{}", i);
    let eos = EOSs::HotQCD;
    let (p, dpde, temp): (Eos, Eos, Eos) = match eos {
        EOSs::WB => (&wb::p, &wb::dpde, &wb::T),
        EOSs::HotQCD => (&hotqcd::p, &hotqcd::dpde, &hotqcd::T),
    };
    println!("{}", name);
    // let init = init_from_entropy_density_2d(t0, s, p, dpde, temp);
    let init = init_from_freestream_2d(t0, s, p, dpde);
    viscous2d::<V, S>(
        &name,
        maxdt,
        er,
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
        freezeout_temp_mev,
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
    r: Scheme<S>,
    nb_trento: usize,
) {
    let trentos = prepare_trento_freestream::<V>(nb_trento);
    let dx = l / V as f64;
    println!("{}", r.name);
    let erpow = 2;
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
        converge(dtmax, dtmin, erpow, |dt, er| {
            hydro2d::<V, S>(
                t0,
                tend,
                dx,
                dt,
                er,
                etaovers,
                zetaovers,
                tempcut,
                freezeout_temp_mev,
                r,
                trento,
            )
        });
    }
}
pub fn run_2d<const V: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_mev: f64,
    nb_trento: usize,
) {
    run_convergence_2d::<V, 1>(
        t0,
        tend,
        l,
        dtmin,
        dtmax,
        etaovers,
        zetaovers,
        tempcut,
        freezeout_temp_mev,
        gauss_legendre_1(Some((0, 1e-3))),
        nb_trento,
    );
    run_convergence_2d::<V, 2>(
        t0,
        tend,
        l,
        dtmin,
        dtmax,
        etaovers,
        zetaovers,
        tempcut,
        freezeout_temp_mev,
        heun(),
        nb_trento,
    );
}
pub fn run_trento_2d<const V: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dt: f64,
    etaovers: (f64, f64, f64),
    zetaovers: (f64, f64, f64),
    tempcut: f64,
    freezeout_temp_mev: f64,
    nb_trento: usize,
) {
    let trentos = prepare_trento_freestream::<V>(nb_trento);
    let gl1 = gauss_legendre_1(Some((0, 1e-3)));
    let heun = heun();
    let dx = l / V as f64;
    let er = dt * dt;
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
        hydro2d::<V, 1>(
            t0,
            tend,
            dx,
            dt,
            er,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_mev,
            gl1,
            trento,
        );
        hydro2d::<V, 2>(
            t0,
            tend,
            dx,
            dt,
            er,
            etaovers,
            zetaovers,
            tempcut,
            freezeout_temp_mev,
            heun,
            trento,
        );
    }
}
