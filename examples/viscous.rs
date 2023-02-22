use std::thread;

use implhydro::{
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
    r: Scheme<S>,
    init_e: ([[f64; V]; V], usize),
) -> HydroOutput<V, V, F_SHEAR_2D, C_SHEAR_2D> {
    let (es, i) = init_e;
    let name = format!("InitTrento{}", i);
    let (p, dpde, temp): (Eos, Eos, Eos) = (&wb::p, &wb::dpde, &wb::T);
    println!("{}", name);
    let init = init_from_energy_2d(t0, es, p, dpde);
    shear2d::<V, S>(
        &name, maxdt, er, t0, tend, dx, r, p, dpde, temp, &init, etaovers,
    )
}

pub fn run_convergence<const V: usize, const S: usize, const TRENTO: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    ermin: f64,
    etaovers: f64,
    r: Scheme<S>,
) {
    let trentos = prepare_trento::<V, TRENTO>();
    let dx = 2.0 * l / V as f64;
    let er0 = (dx * 0.1).powf(2.0); // er0 is so that dt0 = sq2(er0) = dx/2
    let sq2 = |v: f64| v.powf(0.5);
    println!("{}", r.name);
    for i in 0..TRENTO {
        let trento = (trentos[i], i);
        converge(er0, ermin, |er| {
            hydro2d::<V, S>(t0, tend, dx, sq2(er), er, etaovers, r, trento)
        });
    }
}
pub fn run<const V: usize, const TRENTO: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    ermin: f64,
    etaovers: f64,
) {
    run_convergence::<V, 1, TRENTO>(t0, tend, l, ermin, etaovers, gauss_legendre_1(Some(0)));
    run_convergence::<V, 2, TRENTO>(t0, tend, l, ermin, etaovers, heun());
}
pub fn run_trento<const V: usize, const TRENTO: usize>(t0: f64, tend: f64, l: f64, etaovers: f64) {
    let trentos = prepare_trento::<V, TRENTO>();
    // let r = gauss_legendre_1(Some(0));
    let r = heun();
    const S: usize = 2;
    let dx = 2.0 * l / V as f64;
    let dt = dx * 0.1;
    let er = dt * dt;
    for i in 0..TRENTO {
        let trento = (trentos[i], i);
        hydro2d::<V, S>(t0, tend, dx, dt, er, etaovers, r, trento);
    }
}

fn big_stack() {
    let t0 = 1.0;
    let tend = t0 + 3.5;
    let l = 10.0;
    let etaovers = 0.08;

    // let ermin = 1e-5;
    // run::<100, 1>(t0, tend, l, ermin, etaovers);
    // run::<200, 2>(t0, tend, l, ermin, etaovers);

    run_trento::<100, 3>(t0, tend, l, etaovers);
    // run_trento::<200, 100>(t0, tend, l, etaovers);
}

fn main() {
    const STACK_SIZE: usize = 128 * 1024 * 1024; // if you want to run 2D simulation with more than 200x200 cells, you will need to increase the stack size
    thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(big_stack)
        .unwrap()
        .join()
        .unwrap();
}
