use crate::{
    hydro::{
        eos::wb,
        ideal::ideal3d::{self, init_from_energy_density_3d},
        utils::{converge, prepare_trento_3d},
        Eos, HydroOutput, C_IDEAL_3D, F_IDEAL_3D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro3d<const XY: usize, const Z: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    r: Scheme<S>,
    init_e: (&[[[f64; XY]; XY]; Z], usize),
    save_raw: Option<f64>,
) -> HydroOutput<XY, XY, Z, F_IDEAL_3D, C_IDEAL_3D> {
    let (es, i) = init_e;

    let name = ("InitTrento", i);
    let (p, dpde, _temp): (Eos, Eos, Eos) = (&wb::p, &wb::dpde, &wb::T);
    println!("{}{}", name.0, name.1);
    let init = init_from_energy_density_3d(t0, es, p, dpde);
    ideal3d::ideal3d::<XY, Z, S>(&name, maxdt, t0, tend, dx, r, p, dpde, &init, save_raw)
}

pub fn run_convergence_3d<const XY: usize, const Z: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    r: Scheme<S>,
    nb_trento: usize,
    save_raw: Option<f64>,
) {
    let trentos = prepare_trento_3d::<XY, Z>(nb_trento);
    let dx = l / XY as f64;
    let dt0 = dtmax;
    println!("{}", r.name);
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), i);
        converge(dt0, dtmin, |dt| {
            hydro3d::<XY, Z, S>(t0, tend, dx, dt, r, trento, save_raw)
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
    nb_trento: usize,
    save_raw: Option<f64>,
) {
    let imp = gauss_legendre_1();
    const I: usize = 1;
    // let imp = pareschi();
    // const I: usize = 2;
    let exp = heun();
    match solver {
        Solver::Both => {
            run_convergence_3d::<XY, Z, I>(t0, tend, l, dtmin, dtmax, imp, nb_trento, save_raw);
            run_convergence_3d::<XY, Z, 2>(t0, tend, l, dtmin, dtmax, exp, nb_trento, save_raw);
        }
        Solver::Implicit => {
            run_convergence_3d::<XY, Z, I>(t0, tend, l, dtmin, dtmax, imp, nb_trento, save_raw);
        }
        Solver::Explicit => {
            run_convergence_3d::<XY, Z, 2>(t0, tend, l, dtmin, dtmax, exp, nb_trento, save_raw);
        }
    }
}

pub fn run_trento_3d<const XY: usize, const Z: usize>(
    solver: Solver,
    t0: f64,
    tend: f64,
    l: f64,
    dt: f64,
    nb_trento: usize,
    save_raw: Option<f64>,
) {
    let trentos = prepare_trento_3d::<XY, Z>(nb_trento);
    let imp = gauss_legendre_1();
    let exp = heun();
    let dx = l / XY as f64;
    for i in 0..nb_trento {
        let trento = (trentos[i].as_ref(), i);
        match solver {
            Solver::Both => {
                hydro3d::<XY, Z, 1>(t0, tend, dx, dt, imp, trento, save_raw);
                hydro3d::<XY, Z, 2>(t0, tend, dx, dt, exp, trento, save_raw);
            }
            Solver::Implicit => {
                hydro3d::<XY, Z, 1>(t0, tend, dx, dt, imp, trento, save_raw);
            }
            Solver::Explicit => {
                hydro3d::<XY, Z, 2>(t0, tend, dx, dt, exp, trento, save_raw);
            }
        }
    }
}
