use implicit_newton::{
    context::Integration,
    gubser::{gubser_err, init_gubser},
    hydro2d::{hydro2d, Coordinate},
};

pub fn impl2d1<const V: usize>(t0: f64, dx: f64, dt: f64, er: f64) -> ([[[f64; 4]; V]; V], f64) {
    // let r = ([[1.0]], None);
    let r = ([[0.5]], Some([1.0]));
    hydro2d::<V, 1>(
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Coordinate::Milne,
        Integration::FixPoint,
        init_gubser(t0),
    )
}
pub fn impl2d2<const V: usize>(t0: f64, dx: f64, dt: f64, er: f64) -> ([[[f64; 4]; V]; V], f64) {
    // let r = ([[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]], None);
    let sq3 = 3.0f64.sqrt();
    let r = (
        [
            [1.0 / 4.0, 1.0 / 4.0 - 1.0 / 6.0 * sq3],
            [1.0 / 4.0 + 1.0 / 6.0 * sq3, 1.0 / 4.0],
        ],
        Some([1.0 / 2.0, 1.0 / 2.0]),
    );
    hydro2d::<V, 2>(
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Coordinate::Milne,
        Integration::FixPoint,
        init_gubser(t0),
    )
}
pub fn expl2d1<const V: usize>(t0: f64, dx: f64, dt: f64, er: f64) -> ([[[f64; 4]; V]; V], f64) {
    let r = ([[1.0]], None);
    hydro2d::<V, 1>(
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Coordinate::Milne,
        Integration::Explicit,
        init_gubser(t0),
    )
}
pub fn expl2d2<const V: usize>(t0: f64, dx: f64, dt: f64, er: f64) -> ([[[f64; 4]; V]; V], f64) {
    let r = ([[1.0, 0.0], [0.5, 0.5]], None);
    hydro2d::<V, 2>(
        dt,
        er,
        t0,
        t0 + 4.0,
        dx,
        r,
        Coordinate::Milne,
        Integration::Explicit,
        init_gubser(t0),
    )
}

fn main() {
    let t0 = 0.6;
    let dx = 0.1;
    let dt: f64 = 0.0125;
    let er: f64 = (0.1 * dt / (2.0 * dx)).powf(2.0);
    const SIZE: usize = 101;

    // let (v, t) = impl2d1::<SIZE>(t0, dx, dt, er);
    // let (v, t) = impl2d2::<SIZE>(t0, dx, dt, er);
    // let (v, t) = expl2d1::<SIZE>(t0, dx, dt, er);
    let (v, t) = expl2d2::<SIZE>(t0, dx, dt, er);
    let [maxerr, meanerr, maxerrt00, meanerrt00] = gubser_err(v, t, dx);
    println!(
        "|gubser |    e      t00   |\n|-------|-----------------|\n|maxerr | {:.5} {:.5} |\n|meanerr| {:.5} {:.5} |\n",
        maxerr, maxerrt00, meanerr, meanerrt00
    );
}
