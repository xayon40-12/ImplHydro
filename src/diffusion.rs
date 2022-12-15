use crate::{
    newton::{Boundary, Context},
    solver::run,
    utils::{noboundary, periodic},
};

fn flux<const V: usize>(
    v: &[[[f64; 1]; V]; 1],
    boundary: &[Boundary; 2],
    [i, _]: [i32; 2],
) -> [f64; 1] {
    let ivdx = 10.0f64;
    let f = 0;
    let res = ivdx.powf(2.0)
        * (-v[0][boundary[0](i - 2, V)][f] + 16.0 * v[0][boundary[0](i - 1, V)][f]
            - 30.0 * v[0][boundary[0](i, V)][f]
            + 16.0 * v[0][boundary[0](i + 1, V)][f]
            - v[0][boundary[0](i + 2, V)][f])
        / 12.0;
    [res]
}

pub fn diffusion() {
    const V: usize = 100;
    let mut vs = [[[0.0]; V]];
    let k = [[[[0.0]; V]]];
    let integrated = [true];
    for i in 0..V {
        vs[0][i][0] = 100.0 * (-(i as f64 - V as f64 / 2.0).abs()).exp()
    }
    let mean = vs[0].iter().fold(0.0, |acc, i| acc + i[0]) / V as f64;
    let context = Context {
        fun: &flux,
        boundary: &[&periodic, &noboundary], // use noboundary to emulate 1D
        local_interaction: [4, 0],           // use a distance of 0 to emulate 1D
        vs,
        k,
        integrated,
        r: [[1.0]],
        dt: 10.0,
        er: 1e-5,
        tbeg: 0.0,
        tend: 1000.0,
    };
    let (vs, tot_f, tsteps) = run(context);
    let tot_f = tot_f * (2 * 4 + 1 + 1);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
    let mean_res = vs[0].iter().fold(0.0, |acc, i| acc + i[0]) / V as f64;
    println!("expected mean: {}", mean);
    println!("result mean:   {}", mean_res);
}
