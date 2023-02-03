use std::thread;

use implhydro::{
    hydro::{
        eos::wb,
        from_file::{init_from_energy_2d, load_matrix},
        gubser::init_gubser,
        hydro1d, hydro2d, ideal_gas,
        riemann::init_riemann,
        HydroOutput, Pressure, F_IDEAL_1D, F_IDEAL_2D,
    },
    solver::{
        space::Order::{self, *},
        time::schemes::*,
    },
};

pub trait DoHydro<const V: usize, const S: usize, I> {
    type Output;
    fn dohydro(
        name: &str,
        maxdt: f64,
        er: f64,
        t: f64,
        tend: f64,
        dx: f64,
        r: Scheme<S>,
        p: Pressure,
        dpde: Pressure,
        init: I,
        space_order: Order,
    ) -> Self::Output;
}
pub struct Ideal2D;
impl<const V: usize, const S: usize> DoHydro<V, S, Option<([[f64; V]; V], usize)>> for Ideal2D {
    type Output = HydroOutput<V, V, F_IDEAL_2D>;
    fn dohydro(
        name: &str,
        maxdt: f64,
        er: f64,
        t0: f64,
        tend: f64,
        dx: f64,
        r: Scheme<S>,
        p: Pressure,
        dpde: Pressure,
        init_e: Option<([[f64; V]; V], usize)>,
        space_order: Order,
    ) -> Self::Output {
        let init = if let Some((es, _)) = init_e {
            init_from_energy_2d(t0, es, p, dpde)
        } else {
            init_gubser(t0, p, dpde)
        };
        hydro2d::hydro2d::<V, S>(
            name,
            maxdt,
            er,
            t0,
            tend,
            dx,
            r,
            p,
            dpde,
            &init,
            space_order,
        )
    }
}

pub struct Ideal1D;
impl<const V: usize, const S: usize> DoHydro<V, S, bool> for Ideal1D {
    type Output = HydroOutput<V, 1, F_IDEAL_1D>;
    fn dohydro(
        name: &str,
        maxdt: f64,
        er: f64,
        t0: f64,
        tend: f64,
        dx: f64,
        r: Scheme<S>,
        p: Pressure,
        dpde: Pressure,
        use_void: bool,
        space_order: Order,
    ) -> Self::Output {
        hydro1d::hydro1d::<V, S>(
            name,
            maxdt,
            er,
            t0,
            tend,
            dx,
            r,
            p,
            dpde,
            &init_riemann(1.0, p, dpde, use_void),
            space_order,
        )
    }
}

pub fn hydro1d<const V: usize, const S: usize, H: DoHydro<V, S, bool>>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    r: Scheme<S>,
    use_void: bool,
    space_order: Order,
) -> H::Output {
    let void = if use_void { "Void" } else { "" };
    println!("Rieman{}", void);
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    let name = format!("Riemann{}", void);
    H::dohydro(
        &name,
        dt,
        er,
        t0,
        tend,
        dx,
        r,
        &p,
        &dpde,
        use_void,
        space_order,
    )
}

pub fn hydro2d<const V: usize, const S: usize, H: DoHydro<V, S, Option<([[f64; V]; V], usize)>>>(
    t0: f64,
    tend: f64,
    dx: f64,
    dt: f64,
    er: f64,
    r: Scheme<S>,
    init_e: Option<([[f64; V]; V], usize)>,
    space_order: Order,
) -> H::Output {
    let name = if let Some((_, i)) = init_e {
        format!("InitTrento{}", i)
    } else {
        format!("Gubser")
    };
    let (p, dpde): (Pressure, Pressure) = if init_e.is_some() {
        (&wb::p, &wb::dpde)
    } else {
        (&ideal_gas::p, &ideal_gas::dpde)
    };
    println!("{}", name);
    H::dohydro(
        &name,
        dt,
        er,
        t0,
        tend,
        dx,
        r,
        &p,
        &dpde,
        init_e,
        space_order,
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
) -> Option<()> {
    if ermin < 1e-15 {
        eprintln!("ermin<1e-15 in converge, might not converge, set to ermin=1e-15 for safety.");
        ermin = 1e-15;
    }
    let mut f = fun(er)?.0;
    println!("error convergence:");
    er *= 0.1;
    while er > ermin {
        let f2 = fun(er)?.0;
        let (ma, av) = compare(0, &f, &f2);
        println!("er: {:.3e}, max: {:.3e}, average: {:.3e}", er, ma, av);
        f = f2;
        er *= 0.1;
    }
    println!("");

    Some(())
}

pub fn prepare_trento<const V: usize, const TRENTO: usize>() -> [[[f64; V]; V]; TRENTO] {
    let mut trentos = [[[0.0f64; V]; V]; TRENTO];
    for i in 0..TRENTO {
        trentos[i] = load_matrix(&format!("e{}/{:0>2}.dat", V, i)).expect(&format!(
            "Could not load trento initial condition file \"e{}/{:0>2}.dat\".",
            V, i
        ));
    }
    trentos
}

pub fn run<const V: usize, const TRENTO: usize>(
    t0: f64,
    tend: f64,
    dx: f64,
    ermin: f64,
    ord: Order,
) {
    let gl1 = gauss_legendre_1();
    // let gl2 = gauss_legendre_2();
    let heun = heun();
    // let cn = crank_nicolson();

    let trentos = prepare_trento::<V, TRENTO>();

    let er0 = (dx / 2.0).powf(2.0); // er0 is so that dt0 = sq2(er0) = dx/2
    let sq2 = |v: f64| v.powf(0.5);
    // let sq4 = |v: f64| v.powf(0.25);

    let r = gl1;
    println!("{}", r.name);
    converge(er0, ermin, |er| {
        hydro1d::<V, 1, Ideal1D>(t0, tend, dx, sq2(er), er, r, true, ord)
    });
    converge(er0, ermin, |er| {
        hydro1d::<V, 1, Ideal1D>(t0, tend, dx, sq2(er), er, r, false, ord)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 1, Ideal2D>(t0, tend, dx, sq2(er), er, r, None, ord)
    });
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        converge(er0, ermin, |er| {
            hydro2d::<V, 1, Ideal2D>(t0, tend, dx, sq2(er), er, r, trento, ord)
        });
    }
    let r = heun;
    println!("{}", r.name);
    converge(er0, ermin, |er| {
        hydro1d::<V, 2, Ideal1D>(t0, tend, dx, sq2(er), er, r, true, ord)
    });
    converge(er0, ermin, |er| {
        hydro1d::<V, 2, Ideal1D>(t0, tend, dx, sq2(er), er, r, false, ord)
    });
    converge(er0, ermin, |er| {
        hydro2d::<V, 2, Ideal2D>(t0, tend, dx, sq2(er), er, r, None, ord)
    });
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        converge(er0, ermin, |er| {
            hydro2d::<V, 2, Ideal2D>(t0, tend, dx, sq2(er), er, r, trento, ord)
        });
    }
}

fn big_stack() {
    const TRENTO: usize = 100;
    let t0 = 1.0;
    let l = 10.0;
    let tend = t0 + 5.0;
    // let ord = Order2;
    // let ord = Order2Cut(1e-9);
    let ord = Order3;
    // let ord = Order3Cut(1e-9);

    // let ermin = 1e-5;
    // run::<100, TRENTO>(t0, tend, 2.0 * l / 100.0, ermin, ord);
    // run::<200, TRENTO>(t0, tend, 2.0 * l / 200.0, ermin, ord);
    // run::<300, TRENTO>(t0, tend, 2.0 * l / 300.0, ermin, ord);

    const N: usize = 100;
    let trentos = prepare_trento::<N, TRENTO>();
    let gl1 = gauss_legendre_1();
    let dx = 2.0 * l / N as f64;
    let er = 1e-3;
    let sq2 = |v: f64| v.powf(0.5);
    for i in 0..TRENTO {
        let trento = Some((trentos[i], i));
        hydro2d::<N, 1, Ideal2D>(t0, tend, dx, sq2(er), er, gl1, trento, ord);
    }
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
