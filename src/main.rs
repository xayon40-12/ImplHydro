use implicit_newton::{
    diffusion::diffusion,
    hydro1d::hydro1d,
    hydro2d::{hydro2d, Coordinate},
};

fn main() {
    diffusion();
    hydro1d();
    hydro2d::<101>(0.01, 1e-3, 0.6, 4.7, 0.1, Coordinate::Cartesian, |x, y| {
        let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-15 };
        [e, 0.0, 0.0, e]
    });
    // println!("done");
}
