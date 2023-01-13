use std::thread;

use fixhydro::{
    hydro::{
        eos::wb,
        from_file::{init_from_energy_2d, load_matrix},
        gubser::init_gubser,
        hydro1d,
        hydro2d::{self, Coordinate::*},
        ideal_gas,
        riemann::init_riemann,
        Pressure,
    },
    solver::time::schemes::*,
};

pub fn hydro1d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    r: Scheme<S>,
    use_void: bool,
) -> Option<([[[f64; 2]; V]; 1], f64, usize, usize)> {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    hydro1d::hydro1d::<V, S>(
        &format!("Riemann{}", void),
        dt,
        er,
        t0,
        tend,
        dx,
        r,
        p,
        dpde,
        &init_riemann(1.0, p, dpde, use_void),
    )
}
pub fn hydro2d<const V: usize, const S: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    r: Scheme<S>,
    use_exponential: bool,
    init_e: Option<([[f64; V]; V], usize)>,
) -> Option<([[[f64; 3]; V]; V], f64, usize, usize)> {
    let name = if use_exponential {
        if let Some((_, i)) = init_e {
            format!("ExponentialTrento{}", i)
        } else {
            format!("ExponentialGubser")
        }
    } else {
        if let Some((_, i)) = init_e {
            format!("InitTrento{}", i)
        } else {
            format!("Gubser")
        }
    };
    let (p, dpde): (Pressure, Pressure) = if init_e.is_some() {
        (&wb::p, &wb::dpde)
    } else {
        (&ideal_gas::p, &ideal_gas::dpde)
    };
    let init = if let Some((es, _)) = init_e {
        init_from_energy_2d(t0, es, p, dpde)
    } else {
        init_gubser(t0, p, dpde)
    };
    println!("{}", name);
    hydro2d::hydro2d::<V, S>(
        &name,
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
        &init,
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
    mut er: f64,
    mut ermin: f64,
    fun: impl Fn(f64) -> Option<([[[f64; F]; VX]; VY], f64, usize, usize)>,
) -> Option<Vec<[[[f64; F]; VX]; VY]>> {
    if ermin < 1e-15 {
        eprintln!("ermin<1e-15 in converge, might not converge, set to ermin=1e-15 for safety.");
        ermin = 1e-15;
    }
    let mut f = fun(er)?.0;
    let mut all = vec![f];
    println!("error convergence:");
    er *= 0.1;
    while er > ermin {
        let f2 = fun(er)?.0;
        let (ma, av) = compare(0, &f, &f2);
        println!("er: {:.3e}, max: {:.3e}, average: {:.3e}", er, ma, av);
        f = f2;
        all.push(f2);
        er *= 0.1;
    }
    println!("");

    Some(all)
}

pub fn run<const V: usize>(t0: f64, tend: f64, dx: f64, ermin: f64) {
    let gl1 = gauss_legendre_1();
    // let gl2 = gauss_legendre_2();
    let heun = heun();

    const TRENTO: usize = 10;
    let mut trentos = [[[0.0f64; V]; V]; TRENTO];
    for i in 0..TRENTO {
        trentos[i] = load_matrix(&format!("e{}/{:0>2}.dat", V, i)).expect(&format!(
            "Could not load trento initial condition file \"e{}/{:0>2}.dat\".",
            V, i
        ));
    }

    let er0 = (dx / 2.0).powf(2.0); // er0 is so that dt0 = sq2(er0) = dx/2
    let sq2 = |v: f64| v.powf(0.5);
    // let sq4 = |v: f64| v.powf(0.25);

    let r = gl1;
    println!("{}", r.name);
    converge(er0, ermin, |er| {
        hydro1d::<V, 1>(t0, tend, dx, sq2(er), er, r, true)
    });
    converge(er0, ermin, |er| {
        hydro1d::<V, 1>(t0, tend, dx, sq2(er), er, r, false)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 1>(t0, tend, dx, sq2(er), er, r, true, None)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 1>(t0, tend, dx, sq2(er), er, r, false, None)
    });
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        converge(er0, ermin, |er| {
            hydro2d::<V, 1>(t0, tend, dx, sq2(er), er, r, false, trento)
        });
    }
    // let r = gl2;
    // let p = 1.5;
    // println!("{}", r.name);
    // converge(er0.powf(p), ermin.powf(p), |er| {
    //     hydro1d::<V, 2>(t0, tend, dx, sq4(er), er, r, true)
    // });
    // converge(er0.powf(p), ermin.powf(p), |er| {
    //     hydro1d::<V, 2>(t0, tend, dx, sq4(er), er, r, false)
    // });
    // converge(er0.powf(p), ermin.powf(p), |er| {
    //     hydro2d::<V, 2>(t0, tend, dx, sq4(er), er, r, true, None)
    // });
    // converge(er0.powf(p), ermin.powf(p), |er| {
    //     hydro2d::<V, 2>(t0, tend, dx, sq4(er), er, r, false, None)
    // });
    // for i in 0..TRENTO {
    //     let trento = Some((trentos[i], i));
    //     converge(er0.powf(p), ermin.powf(p), |er| {
    //         hydro2d::<V, 2>(t0, tend, dx, sq4(er), er, r, false, trento)
    //     });
    // }
    let r = heun;
    println!("{}", r.name);
    converge(er0, ermin, |er| {
        hydro1d::<V, 2>(t0, tend, dx, sq2(er), er, r, true)
    });
    converge(er0, ermin, |er| {
        hydro1d::<V, 2>(t0, tend, dx, sq2(er), er, r, false)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 2>(t0, tend, dx, sq2(er), er, r, true, None)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 2>(t0, tend, dx, sq2(er), er, r, false, None)
    });
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        converge(er0, ermin, |er| {
            hydro2d::<V, 2>(t0, tend, dx, sq2(er), er, r, false, trento)
        });
    }
}

fn big_stack() {
    let t0 = 1.0;
    let l = 10.0;
    let tend = 4.5;

    let ermin = 1e-9;
    run::<100>(t0, tend, 2.0 * l / 100.0, ermin);
    run::<300>(t0, tend, 2.0 * l / 300.0, ermin);
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
