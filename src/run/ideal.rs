use crate::{
    hydro::{
        eos::wb,
        ideal::{
            ideal1d,
            ideal2d::{self, init_from_entropy_density_2d},
        },
        ideal_gas,
        solutions::{gubser::init_gubser, riemann::init_riemann},
        utils::{converge, prepare_trento},
        Eos, HydroOutput, C_IDEAL_1D, C_IDEAL_2D, F_IDEAL_1D, F_IDEAL_2D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    r: Scheme<S>,
    init_e: Option<([[f64; V]; V], usize)>,
) -> HydroOutput<V, V, F_IDEAL_2D, C_IDEAL_2D> {
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
    ideal2d::ideal2d::<V, S>(&name, maxdt, t0, tend, dx, r, p, dpde, &init)
}

fn hydro1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    r: Scheme<S>,
    use_void: bool,
) -> HydroOutput<V, 1, F_IDEAL_1D, C_IDEAL_1D> {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let name = format!("Riemann{}", void);
    let init = &init_riemann(t0, p, dpde, use_void);
    ideal1d::ideal1d::<V, S>(&name, maxdt, t0, tend, dx, r, p, dpde, &init)
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
    println!("{}", r.name);
    converge(dt0, dtmin, |dt| hydro1d::<V, S>(t0, tend, dx, dt, r, true));
    converge(dt0, dtmin, |dt| hydro1d::<V, S>(t0, tend, dx, dt, r, false));
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
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let dx = l / V as f64;
    let dt0 = dtmax;
    println!("{}", r.name);
    if gubser {
        converge(dt0, dtmin, |dt| hydro2d::<V, S>(t0, tend, dx, dt, r, None));
    }
    for i in 0..nb_trento {
        let trento = Some((trentos[i], i));
        converge(dt0, dtmin, |dt| {
            hydro2d::<V, S>(t0, tend, dx, dt, r, trento)
        });
    }
}
pub fn run_1d<const V: usize>(solver: Solver, t0: f64, tend: f64, l: f64, dtmin: f64, dtmax: f64) {
    // let imp = implicit_euler(); const I: usize = 3;
    // let imp = radauiia2(); const I: usize = 2;
    let imp = gauss_legendre_1();
    const I: usize = 1;
    // let imp = gauss_legendre_2(); const I: usize = 2;
    // let imp = gauss_legendre_3(); const I: usize = 3;
    // let imp = lobatto_iiic();
    // const I: usize = 2;
    // let imp = pareschi();
    // const I: usize = 2;
    let exp = heun();
    const E: usize = 2;
    // let exp = midpoint();
    // const E: usize = 2;
    // let exp = rk4();
    // const E: usize = 4;
    match solver {
        Solver::Both => {
            run_convergence_1d::<V, I>(t0, tend, l, dtmin, dtmax, imp);
            run_convergence_1d::<V, E>(t0, tend, l, dtmin, dtmax, exp);
        }
        Solver::Implicit => {
            run_convergence_1d::<V, I>(t0, tend, l, dtmin, dtmax, imp);
        }
        Solver::Explicit => {
            run_convergence_1d::<V, E>(t0, tend, l, dtmin, dtmax, exp);
        }
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
) {
    let imp = gauss_legendre_1();
    const I: usize = 1;
    // let imp = pareschi();
    // const I: usize = 2;
    let exp = heun();
    match solver {
        Solver::Both => {
            run_convergence_2d::<V, I>(t0, tend, l, dtmin, dtmax, imp, gubser, nb_trento);
            run_convergence_2d::<V, 2>(t0, tend, l, dtmin, dtmax, exp, gubser, nb_trento);
        }
        Solver::Implicit => {
            run_convergence_2d::<V, I>(t0, tend, l, dtmin, dtmax, imp, gubser, nb_trento);
        }
        Solver::Explicit => {
            run_convergence_2d::<V, 2>(t0, tend, l, dtmin, dtmax, exp, gubser, nb_trento);
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
) {
    let trentos = prepare_trento::<V>(nb_trento);
    let imp = gauss_legendre_1();
    let exp = heun();
    let dx = l / V as f64;
    for i in 0..nb_trento {
        let trento = Some((trentos[i], i));
        match solver {
            Solver::Both => {
                hydro2d::<V, 1>(t0, tend, dx, dt, imp, trento);
                hydro2d::<V, 2>(t0, tend, dx, dt, exp, trento);
            }
            Solver::Implicit => {
                hydro2d::<V, 1>(t0, tend, dx, dt, imp, trento);
            }
            Solver::Explicit => {
                hydro2d::<V, 2>(t0, tend, dx, dt, exp, trento);
            }
        }
    }
}
