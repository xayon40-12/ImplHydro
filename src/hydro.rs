use crate::{
    newton::{Boundary, Context},
    solver::run,
    utils::periodic,
};

fn flux<const V: usize>(v: &[[f64; 2]; V], boundary: Boundary, i: i32) -> [f64; 2] {
    let ivdx = 10.0f64;
    let diff = ivdx.powf(2.0)
        * (-v[boundary(i - 2, V)][0] + 16.0 * v[boundary(i - 1, V)][0]
            - 30.0 * v[boundary(i, V)][0]
            + 16.0 * v[boundary(i + 1, V)][0]
            - v[boundary(i + 2, V)][0])
        / 12.0;
    let vi = v[boundary(i, V)][0];
    let svi = v[boundary(i, V)][1];
    let sq = svi + svi * svi - vi; // compute sqrt to check if non integrated function works
    [diff, sq]
}

pub fn hydro() {
    const V: usize = 100;
    let mut vs = [[0.0; V], [0.0; V]];
    let integrated = [true, false];
    for i in 0..V {
        if i < V / 2 {
            vs[0][i] = 10.0;
        } else {
            vs[0][i] = 1.0;
        }
    }
    let mean = vs[0].iter().fold(0.0, |acc, i| acc + i) / V as f64;
    let context = Context {
        fun: &flux,
        boundary: &periodic,
        local_interaction: 4,
        vs,
        integrated,
        r: [[1.0]],
        dt: 1.0,
        er: 1e-10,
        tbeg: 0.0,
        tend: 100.0,
    };
    let (vs, tot_f, tsteps) = run(context);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
    let mean_res = vs[0].iter().fold(0.0, |acc, i| acc + i) / V as f64;
    let mean_res1 = vs[1].iter().fold(0.0, |acc, i| acc + i) / V as f64;
    println!("expected mean: {}", mean);
    println!("result mean:   {}", mean_res);
    println!("expected mean sqrt: {}", mean.sqrt());
    println!("result mean sqrt:   {}", mean_res1);
}
