use implicit_newton::{
    context::Integration::*,
    gubser::{gubser_err, init_gubser},
    hydro1d::hydro1d,
    hydro2d::{hydro2d, Coordinate::Milne},
    riemann::init_riemann,
    schemes::*,
};

pub fn impl1d1<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 3]; V]; 1], f64, usize, usize) {
    let r = gauss_legendre_1();
    hydro1d::<V, 1>(
        "impl1d1",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        FixPoint,
        init_riemann(),
    )
}
pub fn expl1d2<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 3]; V]; 1], f64, usize, usize) {
    let r = heun();
    hydro1d::<V, 2>(
        "expl1d2",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Explicit,
        init_riemann(),
    )
}
pub fn impl2d1<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    // let r = ([[1.0]], None);
    let r = gauss_legendre_1();
    hydro2d::<V, 1>(
        "impl2d1",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Milne,
        FixPoint,
        init_gubser(t0),
    )
}
pub fn impl2d2<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    let r = gauss_legendre_2();
    hydro2d::<V, 2>(
        "impl2d2",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Milne,
        FixPoint,
        init_gubser(t0),
    )
}
pub fn expl2d1<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    let r = euler();
    hydro2d::<V, 1>(
        "expl2d1",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Milne,
        Explicit,
        init_gubser(t0),
    )
}
pub fn expl2d2<const V: usize>(
    t0: f64,
    dx: f64,
    dt: f64,
    er: f64,
) -> ([[[f64; 4]; V]; V], f64, usize, usize) {
    let r = heun();
    hydro2d::<V, 2>(
        "expl2d2",
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Milne,
        Explicit,
        init_gubser(t0),
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
    m: i32,
    fun: impl Fn(f64) -> [[[f64; F]; VX]; VY],
) -> Vec<[[[f64; F]; VX]; VY]> {
    let mut dt = 0.05;
    let mut f = fun(dt);
    let mut all = vec![f];
    for i in 0..=m {
        dt *= 0.5;
        let f2 = fun(dt);
        let (ma, av) = compare(0, &f, &f2);
        let d = 0.5f64.powi(i);
        println!(
            "i: {:0>2}, dt: {:.3e}, d: {:.3e}, ma: {:.3e}, av: {:.3e}, ma/d: {:.3e}, av/d: {:.3e}",
            i,
            dt,
            d,
            ma,
            av,
            ma / d,
            av / d
        );
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

    // let (_v, _t, _maxerr, _meanerr) = impl1d1::<SIZE>(t0, dx, dt, er);
    // let (_v, _t, _maxerr, _meanerr) = expl1d2::<SIZE>(t0, dx, dt, er);
    let _v = converge(8, |dt| expl1d2::<SIZE>(t0, dx, dt, dt * dt).0);
    let _v = converge(6, |dt| expl2d2::<SIZE>(t0, dx, dt, dt * dt).0);
    // let (v, t, _maxerr, _meanerr) = impl2d1::<SIZE>(t0, dx, dt, er);
    // let (v, t, _maxerr, _meanerr) = impl2d2::<SIZE>(t0, dx, dt, er);
    // let (v, t, _maxerr, _meanerr) = expl2d1::<SIZE>(t0, dx, dt, er);
    let (v, t, _maxerr, _meanerr) = expl2d2::<SIZE>(t0, dx, dt, er);
    let [maxerr, meanerr, maxerrt00, meanerrt00] = gubser_err(v, t, dx);
    println!(
        "|gubser |    e      t00   |\n|-------|-----------------|\n|maxerr | {:.5} {:.5} |\n|meanerr| {:.5} {:.5} |\n",
        maxerr, maxerrt00, meanerr, meanerrt00
    );
}
