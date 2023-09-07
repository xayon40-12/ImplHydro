use crate::{
    hydro::{
        eos::wb,
        ideal::ideal2d::{self, init_from_entropy_density_2d},
        ideal_gas,
        solutions::gubser::init_gubser,
        utils::{converge, prepare_trento_2d},
        Eos, HydroOutput, C_IDEAL_2D, F_IDEAL_2D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    r: Scheme<S>,
    init_e: Option<(&[[f64; V]; V], usize)>,
    save_raw: bool,
) -> HydroOutput<V, V, 1, F_IDEAL_2D, C_IDEAL_2D> {
    let name = if let Some((_, i)) = init_e {
        format!("InitTrento{}", i)
    } else {
        format!("Gubser")
    };
    let (p, dpde, _temp): (Eos, Eos, Eos) = if init_e.is_some() {
        (&wb::p, &wb::dpde, &wb::T)
    } else {
        (&ideal_gas::p, &ideal_gas::dpde, &ideal_gas::T)
    };
    println!("{}", name);
    let init = if let Some((es, _)) = init_e {
        init_from_entropy_density_2d(t0, es, p, dpde)
    } else {
        init_gubser(t0, p, dpde)
    };
    ideal2d::ideal2d::<V, S>(&name, maxdt, t0, tend, dx, r, p, dpde, &init, save_raw)
}

pub fn run_convergence_2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    r: Scheme<S>,
    gubser: bool,
    nb_trento: usize,
    save_raw: bool,
) {
    let trentos = prepare_trento_2d::<V>(nb_trento);
    let dx = l / V as f64;
    let dt0 = dtmax;
    println!("{}", r.name);
    if gubser {
        converge(dt0, dtmin, |dt| {
            hydro2d::<V, S>(t0, tend, dx, dt, r, None, save_raw)
        });
    }
    for i in 0..nb_trento {
        let trento = Some((trentos[i].as_ref(), i));
        converge(dt0, dtmin, |dt| {
            hydro2d::<V, S>(t0, tend, dx, dt, r, trento, save_raw)
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
    gubser: bool,
    nb_trento: usize,
    save_raw: bool,
) {
    let imp = gauss_legendre_1();
    const I: usize = 1;
    // let imp = pareschi();
    // const I: usize = 2;
    let exp = heun();
    match solver {
        Solver::Both => {
            run_convergence_2d::<V, I>(t0, tend, l, dtmin, dtmax, imp, gubser, nb_trento, save_raw);
            run_convergence_2d::<V, 2>(t0, tend, l, dtmin, dtmax, exp, gubser, nb_trento, save_raw);
        }
        Solver::Implicit => {
            run_convergence_2d::<V, I>(t0, tend, l, dtmin, dtmax, imp, gubser, nb_trento, save_raw);
        }
        Solver::Explicit => {
            run_convergence_2d::<V, 2>(t0, tend, l, dtmin, dtmax, exp, gubser, nb_trento, save_raw);
        }
    }
}

pub fn run_trento_2d<const V: usize>(
    solver: Solver,
    t0: f64,
    tend: f64,
    l: f64,
    dt: f64,
    nb_trento: usize,
    save_raw: bool,
) {
    let trentos = prepare_trento_2d::<V>(nb_trento);
    let imp = gauss_legendre_1();
    let exp = heun();
    let dx = l / V as f64;
    for i in 0..nb_trento {
        let trento = Some((trentos[i].as_ref(), i));
        match solver {
            Solver::Both => {
                hydro2d::<V, 1>(t0, tend, dx, dt, imp, trento, save_raw);
                hydro2d::<V, 2>(t0, tend, dx, dt, exp, trento, save_raw);
            }
            Solver::Implicit => {
                hydro2d::<V, 1>(t0, tend, dx, dt, imp, trento, save_raw);
            }
            Solver::Explicit => {
                hydro2d::<V, 2>(t0, tend, dx, dt, exp, trento, save_raw);
            }
        }
    }
}
