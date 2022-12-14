use crate::{
    newton::{Boundary, Fun},
    solver::run,
    utils::periodic,
};

fn flux<const F: usize, const V: usize>(v: &[[f64; F]; V], boundary: Boundary, i: i32) -> [f64; F] {
    let ivdx = 10.0f64;
    let mut res = [0.0; F];
    for f in 0..F {
        res[f] = ivdx.powf(2.0)
            * (-v[boundary(i - 2, V)][f] + 16.0 * v[boundary(i - 1, V)][f]
                - 30.0 * v[boundary(i, V)][f]
                + 16.0 * v[boundary(i + 1, V)][f]
                - v[boundary(i + 2, V)][f])
            / 12.0
    }
    res
}

pub fn diffusion() {
    const V: usize = 100;
    const F: usize = 1;
    let f: (Fun<F, V>, Boundary, i32) = (&flux, &periodic, 4);
    let mut vs = [[0.0; V]; F];
    for i in 0..V {
        vs[0][i] = 100.0 * (-((i - V / 2) as f64).abs()).exp()
    }
    let tot = vs[0].iter().fold(0.0, |acc, i| acc + i);
    //let r = vec![vec![5.0 / 12.0, -1.0 / 12.0], vec![3.0 / 4.0, 1.0 / 4.0]];
    let r = [[1.0]];
    let dt = 0.1;
    let tbeg = 0.0;
    let tend = 10.0;
    run(f, vs, r, dt, tbeg, tend);
    println!("expected mean: {}", tot / V as f64);
}
