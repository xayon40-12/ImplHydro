use crate::{
    hydro::{
        ideal::ideal1d, ideal_gas, solutions::riemann::init_riemann, utils::converge, HydroOutput,
        C_IDEAL_1D, F_IDEAL_1D,
    },
    solver::{time::schemes::*, Solver},
};

fn hydro1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    maxdt: f64,
    r: Scheme<S>,
    use_void: bool,
    save_raw: Option<f64>,
) -> HydroOutput<V, 1, 1, F_IDEAL_1D, C_IDEAL_1D> {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let name = format!("Riemann{}", void);
    let init = &init_riemann(t0, p, dpde, use_void);
    ideal1d::ideal1d::<V, S>(
        &(&name, 0),
        maxdt,
        t0,
        tend,
        dx,
        r,
        p,
        dpde,
        &init,
        save_raw,
    )
}
pub fn run_convergence_1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    r: Scheme<S>,
    save_raw: Option<f64>,
) {
    let dx = l / V as f64;
    let dt0 = dtmax;
    println!("{}", r.name);
    converge(dt0, dtmin, |dt| {
        hydro1d::<V, S>(t0, tend, dx, dt, r, true, save_raw)
    });
    converge(dt0, dtmin, |dt| {
        hydro1d::<V, S>(t0, tend, dx, dt, r, false, save_raw)
    });
}
pub fn run_1d<const V: usize>(
    solver: Solver,
    t0: f64,
    tend: f64,
    l: f64,
    dtmin: f64,
    dtmax: f64,
    save_raw: Option<f64>,
) {
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
            run_convergence_1d::<V, I>(t0, tend, l, dtmin, dtmax, imp, save_raw);
            run_convergence_1d::<V, E>(t0, tend, l, dtmin, dtmax, exp, save_raw);
        }
        Solver::Implicit => {
            run_convergence_1d::<V, I>(t0, tend, l, dtmin, dtmax, imp, save_raw);
        }
        Solver::Explicit => {
            run_convergence_1d::<V, E>(t0, tend, l, dtmin, dtmax, exp, save_raw);
        }
    }
}
