use crate::{
    hydro::{
        eos::wb,
        utils::{converge, prepare_trento},
        viscous::init_from_energy_2d,
        viscous::shear2d::shear2d,
        Eos, HydroOutput, C_SHEAR_2D, F_SHEAR_2D,
    },
    solver::time::schemes::*,
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    er: f64,
    etaovers: f64,
    tempcut: f64,
    r: Scheme<S>,
    init_e: ([[f64; V]; V], usize),
) -> HydroOutput<V, V, F_SHEAR_2D, C_SHEAR_2D> {
    let (es, i) = init_e;
    let name = format!("InitTrento{}", i);
    let (p, dpde, temp): (Eos, Eos, Eos) = (&wb::p, &wb::dpde, &wb::T);
    println!("{}", name);
    let init = init_from_energy_2d(t0, es, p, dpde);
    shear2d::<V, S>(
        &name, maxdt, er, t0, tend, dx, r, p, dpde, temp, &init, etaovers, tempcut,
    )
}

pub fn run_convergence<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    etaovers: f64,
    tempcut: f64,
    r: Scheme<S>,
    nb_trento: usize,
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let dx = l / V as f64;
    println!("{}", r.name);
    let erpow = 2;
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
        converge(dtmax, dtmin, erpow, |dt, er| {
            hydro2d::<V, S>(t0, tend, dx, dt, er, etaovers, tempcut, r, trento)
        });
    }
}
pub fn run<const V: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    etaovers: f64,
    tempcut: f64,
    nb_trento: usize,
) {
    run_convergence::<V, 1>(
        t0,
        tend,
        l,
        dtmin,
        dtmax,
        etaovers,
        tempcut,
        gauss_legendre_1(Some(0)),
        nb_trento,
    );
    run_convergence::<V, 2>(
        t0,
        tend,
        l,
        dtmin,
        dtmax,
        etaovers,
        tempcut,
        heun(),
        nb_trento,
    );
}
pub fn run_trento<const V: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dt: f64,
    etaovers: f64,
    tempcut: f64,
    nb_trento: usize,
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let r = gauss_legendre_1(Some(0));
    // let r = heun();
    const S: usize = 1;
    let dx = l / V as f64;
    let er = dt * dt;
    for i in 0..nb_trento {
        let trento = (trentos[i], i);
        hydro2d::<V, S>(t0, tend, dx, dt, er, etaovers, tempcut, r, trento);
    }
}
