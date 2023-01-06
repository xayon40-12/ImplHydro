use std::thread;

use fixhydro::{
    hydro::{
        gubser::init_gubser,
        hydro1d,
        hydro2d::{self, Coordinate::*},
        ideal_gas,
        riemann::init_riemann,
        Pressure,
    },
    solver::schemes::*,
};

pub fn hydro1d<'a, const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
    r: Scheme<S>,
    use_void: bool,
) -> ([[[f64; 3]; V]; 1], f64, usize, usize) {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    hydro1d::hydro1d::<V, S>(
        &format!("Riemann{}", void),
        dt,
        er,
        t0,
        tend,
        dx,
        r,
        &p,
        &dpde,
        init_riemann(&p, &dpde, use_void),
    )
}
pub fn hydro2d<'a, const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
    r: Scheme<S>,
    use_exponential: bool,
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    let name = if use_exponential {
        "Exponential"
    } else {
        "Gubser"
    };
    println!("{}", name);
    hydro2d::hydro2d::<V, S>(
        name,
        dt,
        er,
        t0,
        tend,
        dx,
        r,
        Milne,
        // Cartesian,
        &p,
        &dpde,
        init_gubser(t0, &p, &dpde),
        use_exponential,
    )
}

pub fn compare<const VX: usize, const VY: usize, const F: usize>(
    f: usize,
    v1: &[[[f64; F]; VX]; VY],
    v2: &[[[f64; F]; VX]; VY],
) -> (f64, f64) {
    let mut average: f64 = 0.0;
    let mut maximum: f64 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            let v1 = v1[j][i][f];
            let v2 = v2[j][i][f];
            let d = (v1 - v2).abs() / v1.abs().max(v2.abs());
            average += d;
            maximum = maximum.max(d);
        }
    }
    average /= (VX * VY) as f64;

    (maximum, average)
}

pub fn converge<const VX: usize, const VY: usize, const F: usize>(
    mut dt: f64,
    m: usize,
    fun: impl Fn(f64) -> [[[f64; F]; VX]; VY],
) -> Vec<[[[f64; F]; VX]; VY]> {
    let mut f = fun(dt);
    let mut all = vec![f];
    println!("error convergence:");
    for _i in 0..m {
        dt *= 0.5;
        let f2 = fun(dt);
        let (ma, av) = compare(0, &f, &f2);
        println!("dt: {:.3e}, max: {:.3e}, average: {:.3e}", dt, ma, av);
        f = f2;
        all.push(f2);
    }
    println!("");

    all
}

pub fn run<const V: usize>(t0: f64, tend: f64, dx: f64, nconv: usize) {
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let gl1 = gauss_legendre_1();
    let gl2 = gauss_legendre_2();
    let heun = heun();

    // let r = gl1;
    // let dt = dx;
    // println!("{}", r.name);
    // converge(dt, nconv, |dt| {
    //     hydro1d::<V, 1>(t0, tend, dx, dt, dt * dt, p, dpde, r, true).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro1d::<V, 1>(t0, tend, dx, dt, dt * dt, p, dpde, r, false).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro2d::<V, 1>(t0, tend, dx, dt, dt * dt, p, dpde, r, true).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro2d::<V, 1>(t0, tend, dx, dt, dt * dt, p, dpde, r, false).0
    // });
    // let r = gl2;
    // let dt = dx * 4.0;
    // println!("{}", r.name);
    // converge(dt, nconv, |dt| {
    //     hydro1d::<V, 2>(t0, tend, dx, dt, dt * dt * dt * dt, p, dpde, r, true).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro1d::<V, 2>(t0, tend, dx, dt, dt * dt * dt * dt, p, dpde, r, false).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro2d::<V, 2>(t0, tend, dx, dt, dt * dt * dt * dt, p, dpde, r, true).0
    // });
    // converge(dt, nconv, |dt| {
    //     hydro2d::<V, 2>(t0, tend, dx, dt, dt * dt * dt * dt, p, dpde, r, false).0
    // });
    let r = heun;
    let dt = dx / 2.0;
    println!("{}", r.name);
    converge(dt, nconv, |dt| {
        println!("dt: {}", dt);
        hydro1d::<V, 2>(t0, tend, dx, dt, dt * dt, p, dpde, r, true).0
    });
    converge(dt, nconv, |dt| {
        hydro1d::<V, 2>(t0, tend, dx, dt, dt * dt, p, dpde, r, false).0
    });
    converge(dt, nconv, |dt| {
        hydro2d::<V, 2>(t0, tend, dx, dt, dt * dt, p, dpde, r, true).0
    });
    converge(dt, nconv, |dt| {
        hydro2d::<V, 2>(t0, tend, dx, dt, dt * dt, p, dpde, r, false).0
    });
}

fn big_stack() {
    let t0 = 1.0;
    let nconv = 2;

    let tend = 4.5;
    run::<100>(t0, tend, 0.1, nconv);
    // run::<200>(t0, tend, 0.05, nconv);
    // run::<500>(t0, tend, 0.01, nconv);

    // let tend = 9.0;
    // run::<100>(t0, tend, 0.2, nconv);
    // run::<200>(t0, tend, 0.1, nconv);
}

fn main() {
    const STACK_SIZE: usize = 16 * 1024 * 1024; // if you want to run 2D simulation with more than 200x200 cells, you will need to increase the stack size
    thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(big_stack)
        .unwrap()
        .join()
        .unwrap();
}
