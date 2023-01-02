use implicit_newton::{
    hydro::{
        gubser::{gubser_err, init_gubser},
        hydro1d,
        hydro2d::{self, Coordinate::Milne},
        ideal_gas,
        riemann::init_riemann,
        Pressure,
    },
    solver::{
        context::Integration::{self, *},
        schemes::*,
    },
};

pub fn hydro1d<'a, const V: usize, const S: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn(Integration, Scheme<S>) -> ([[[f64; 3]; V]; 1], f64, usize, usize) + 'a> {
    Box::new(move |explimpl, r| {
        hydro1d::hydro1d::<V, S>(
            &format!("{:?}1d{}_{:e}", explimpl, S, dt),
            dt,
            er,
            t0,
            t0 + 4.0,
            dx,
            r,
            explimpl,
            &p,
            &dpde,
            init_riemann(&p, &dpde),
        )
    })
}
pub fn hydro2d<'a, const V: usize, const S: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn(Integration, Scheme<S>) -> ([[[f64; 4]; V]; V], f64, usize, usize) + 'a> {
    Box::new(move |explimpl, r| {
        hydro2d::hydro2d::<V, S>(
            &format!("{:?}2d{}_{:e}", explimpl, S, dt),
            dt,
            er,
            t0,
            t0 + 4.0,
            dx,
            r,
            Milne,
            explimpl,
            &p,
            &dpde,
            init_gubser(t0, &p, &dpde),
        )
    })
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
    m: i32,
    fun: impl Fn(f64) -> [[[f64; F]; VX]; VY],
) -> Vec<[[[f64; F]; VX]; VY]> {
    let mut dt = 0.05;
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

fn main() {
    let t0 = 0.6;
    let dx = 0.1;
    let dt: f64 = 0.0125;
    let er: f64 = dt * dt;
    const SIZE: usize = 101;
    let p = &ideal_gas::p;
    let dpde = &ideal_gas::dpde;
    // let hydro1d1 = hydro1d::<SIZE, 1>(t0, dx, dt, er, p, dpde);
    // let hydro1d2 = hydro1d::<SIZE, 2>(t0, dx, dt, er, p, dpde);
    let hydro2d1 = hydro2d::<SIZE, 1>(t0, dx, dt, er, p, dpde);
    // let hydro2d2 = hydro2d::<SIZE, 2>(t0, dx, dt, er, p, dpde);

    // Convergenc:
    // let _v = converge(8, |dt| {
    //     hydro1d::<SIZE, 2>(t0, dx, dt, er, p, dpde)(Explicit, heun()).0
    // });
    // let _v = converge(8, |dt| {
    //     hydro1d::<SIZE, 1>(t0, dx, dt, er, p, dpde)(FixPoint, gauss_legendre_1()).0
    // });

    // let _v = converge(6, |dt| {
    //     hydro2d::<SIZE, 2>(t0, dx, dt, er, p, dpde)(Explicit, heun()).0
    // });
    // let _v = converge(6, |dt| {
    //     hydro2d::<SIZE, 1>(t0, dx, dt, er, p, dpde)(FixPoint, gauss_legendre_1()).0
    // });

    // let (v, t, cost, tsteps) = hydro2d1(Explicit, euler());
    // let (v, t, cost, tsteps) = hydro2d1(FixPoint, euler());
    let (v, t, cost, tsteps) = hydro2d1(FixPoint, gauss_legendre_1());
    // let (v, t, cost, tsteps) = hydro2d2(Explicit, heun());
    // let (v, t, cost, tsteps) = hydro2d2(FixPoint, gauss_legendre_2());
    let [maxerrt00, meanerrt00] = gubser_err(v, t, dx, p);
    println!("cost: {}, tsteps: {}", cost, tsteps);
    println!(
        "|gubser |   t00   |\n|-------|---------|\n|maxerr | {:.5} |\n|meanerr| {:.5} |\n",
        maxerrt00, meanerrt00
    );
}
