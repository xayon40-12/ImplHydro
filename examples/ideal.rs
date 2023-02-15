use std::thread;

use implhydro::{
    hydro::{
        eos::wb,
        ideal::init_from_energy_2d,
        ideal::{ideal1d, ideal2d},
        ideal_gas,
        solutions::{gubser::init_gubser, riemann::init_riemann},
        utils::{converge, prepare_trento},
        Eos, HydroOutput, F_IDEAL_1D, F_IDEAL_2D,
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
) -> HydroOutput<V, V, F_IDEAL_2D> {
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
) -> HydroOutput<V, 1, F_IDEAL_1D> {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let name = format!("Riemann{}", void);
    let init = &init_riemann(1.0, p, dpde, use_void);
    ideal1d::ideal1d::<V, S>(&name, maxdt, er, t0, tend, dx, r, p, dpde, &init)
}

pub fn run_convergence<const V: usize, const S: usize, const TRENTO: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    ermin: f64,
    r: Scheme<S>,
) {
    let trentos = prepare_trento::<V, TRENTO>();
    let dx = 2.0 * l / V as f64;
    let er0 = (dx / 2.0).powf(2.0); // er0 is so that dt0 = sq2(er0) = dx/2
    let sq2 = |v: f64| v.powf(0.5);
    println!("{}", r.name);
    converge(er0, ermin, |er| {
        hydro1d::<V, S>(t0, tend, dx, sq2(er), er, r, true)
    });
    converge(er0, ermin, |er| {
        hydro1d::<V, S>(t0, tend, dx, sq2(er), er, r, false)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, S>(t0, tend, dx, sq2(er), er, r, None)
    });
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        converge(er0, ermin, |er| {
            hydro2d::<V, S>(t0, tend, dx, sq2(er), er, r, trento)
        });
    }
}
pub fn run<const V: usize, const TRENTO: usize>(t0: f64, tend: f64, l: f64, ermin: f64) {
    run_convergence::<V, 1, TRENTO>(t0, tend, l, ermin, gauss_legendre_1());
    run_convergence::<V, 2, TRENTO>(t0, tend, l, ermin, heun());
}
pub fn run_trento<const V: usize, const TRENTO: usize>(t0: f64, tend: f64, l: f64) {
    let trentos = prepare_trento::<V, TRENTO>();
    let gl1 = gauss_legendre_1();
    let dx = 2.0 * l / V as f64;
    let dt = dx * 0.1;
    let er = dt * dt;
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        hydro2d::<V, 1>(t0, tend, dx, dt, er, gl1, trento);
    }
}

fn big_stack() {
    let t0 = 1.0;
    let tend = t0 + 3.5;
    let l = 10.0;

    let ermin = 1e-4;
    run::<100, 2>(t0, tend, l, ermin);
    run::<200, 2>(t0, tend, l, ermin);

    run_trento::<100, 100>(t0, tend, l);
    run_trento::<200, 100>(t0, tend, l);
}

fn main() {
    const STACK_SIZE: usize = 64 * 1024 * 1024; // if you want to run 2D simulation with more than 200x200 cells, you will need to increase the stack size
    thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(big_stack)
        .unwrap()
        .join()
        .unwrap();
}
