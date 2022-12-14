use crate::{
    newton::{Boundary, Context},
    solver::run,
    utils::periodic,
};

fn flux<const V: usize>(v: &[[f64; 1]; V], boundary: Boundary, i: i32) -> [f64; 1] {
    let ivdx = 10.0f64;
    let f = 0;
    let res = ivdx.powf(2.0)
        * (-v[boundary(i - 2, V)][f] + 16.0 * v[boundary(i - 1, V)][f]
            - 30.0 * v[boundary(i, V)][f]
            + 16.0 * v[boundary(i + 1, V)][f]
            - v[boundary(i + 2, V)][f])
        / 12.0;
    [res]
}

pub fn diffusion() {
    const V: usize = 100;
    let mut vs = [[0.0; V]];
    let integrated = [true];
    for i in 0..V {
        vs[0][i] = 100.0 * (-((i - V / 2) as f64).abs()).exp()
    }
    let mean = vs[0].iter().fold(0.0, |acc, i| acc + i) / V as f64;
    let context = Context {
        fun: &flux,
        boundary: &periodic,
        local_interaction: 4,
        vs,
        integrated,
        r: [[1.0]],
        dt: 10.0,
        er: 1e-10,
        tbeg: 0.0,
        tend: 1000.0,
    };
    let (vs, tot_f, tsteps) = run(context);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
    let mean_res = vs[0].iter().fold(0.0, |acc, i| acc + i) / V as f64;
    println!("expected mean: {}", mean);
    println!("result mean:   {}", mean_res);
}
