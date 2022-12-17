use implicit_newton::{
    diffusion::diffusion,
    hydro1d::hydro1d,
    hydro2d::{hydro2d, Coordinate},
};

fn main() {
    // diffusion();
    // hydro1d();
    hydro2d(0.01, 1e-5, 0.6, 4.6, Coordinate::Milne);
    // println!("done");
}
