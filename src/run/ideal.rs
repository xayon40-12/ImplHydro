use crate::{
    hydro::{
        eos::wb,
        ideal::init_from_energy_2d,
        ideal::{ideal1d, ideal2d},
        ideal_gas,
        solutions::{gubser::init_gubser, riemann::init_riemann},
        utils::{converge, prepare_trento},
        Eos, HydroOutput, C_IDEAL_1D, C_IDEAL_2D, F_IDEAL_1D, F_IDEAL_2D,
    },
    solver::time::schemes::*,
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    er: f64,
    r: Scheme<S>,
    init_e: Option<([[f64; V]; V], usize)>,
) -> HydroOutput<V, V, F_IDEAL_2D, C_IDEAL_2D> {
    let name = if let Some((_, i)) = init_e {
        format!("InitTrento{}", i)
    } else {
        format!("Gubser")
    };
    let (p, dpde): (Eos, Eos) = if init_e.is_some() {
        (&wb::p, &wb::dpde)
    } else {
        (&ideal_gas::p, &ideal_gas::dpde)
    };
    println!("{}", name);
    let init = if let Some((es, _)) = init_e {
        init_from_energy_2d(t0, es, p, dpde)
    } else {
        init_gubser(t0, p, dpde)
    };
    ideal2d::ideal2d::<V, S>(&name, maxdt, er, t0, tend, dx, r, p, dpde, &init)
}

fn hydro1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    er: f64,
    r: Scheme<S>,
    use_void: bool,
) -> HydroOutput<V, 1, F_IDEAL_1D, C_IDEAL_1D> {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let name = format!("Riemann{}", void);
    let init = &init_riemann(t0, p, dpde, use_void);
    ideal1d::ideal1d::<V, S>(&name, maxdt, er, t0, tend, dx, r, p, dpde, &init)
}

pub fn run_convergence_1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    r: Scheme<S>,
) {
    let dx = l / V as f64;
    let dt0 = dtmax;
    let erpow = 2;
    println!("{}", r.name);
    converge(dt0, dtmin, erpow, |dt, er| {
        hydro1d::<V, S>(t0, tend, dx, dt, er, r, true)
    });
    converge(dt0, dtmin, erpow, |dt, er| {
        hydro1d::<V, S>(t0, tend, dx, dt, er, r, false)
    });
}

pub fn run_convergence_2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    r: Scheme<S>,
    nb_trento: usize,
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let dx = l / V as f64;
    let dt0 = dtmax;
    let erpow = 2;
    println!("{}", r.name);
    converge(dt0, dtmin, erpow, |dt, er| {
        hydro2d::<V, S>(t0, tend, dx, dt, er, r, None)
    });
    for i in 0..nb_trento {
        let trento = Some((trentos[i], i));
        converge(dt0, dtmin, erpow, |dt, er| {
            hydro2d::<V, S>(t0, tend, dx, dt, er, r, trento)
        });
    }
}
pub fn run_1d<const V: usize>(t0: f64, tend: f64, l: f64, dtmin: f64, dtmax: f64) {
    run_convergence_1d::<V, 1>(t0, tend, l, dtmin, dtmax, gauss_legendre_1(Some((0, 1e-2))));
    run_convergence_1d::<V, 2>(t0, tend, l, dtmin, dtmax, heun());
}
pub fn run_2d<const V: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    nb_trento: usize,
) {
    run_convergence_2d::<V, 1>(
        t0,
        tend,
        l,
        dtmin,
        dtmax,
        gauss_legendre_1(Some((0, 1e-2))),
        nb_trento,
    );
    run_convergence_2d::<V, 2>(t0, tend, l, dtmin, dtmax, heun(), nb_trento);
}
pub fn run_trento_2d<const V: usize>(t0: f64, tend: f64, l: f64, dt: f64, nb_trento: usize) {
    let trentos = prepare_trento::<V>(nb_trento);
    let gl1 = gauss_legendre_1(Some((0, 1e-2)));
    let heun = heun();
    let dx = l / V as f64;
    let er = dt * dt;
    for i in 0..nb_trento {
        let trento = Some((trentos[i], i));
        hydro2d::<V, 1>(t0, tend, dx, dt, er, gl1, trento);
        hydro2d::<V, 2>(t0, tend, dx, dt, er, heun, trento);
    }
}
