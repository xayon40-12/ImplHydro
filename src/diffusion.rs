use crate::{
    context::{Boundary, Context, Integration, ToCompute},
    solver::run,
    utils::{noboundary, periodic},
};

fn flux<const V: usize>(
    [_vo, v]: [&[[[f64; 1]; V]; 1]; 2],
    boundary: &[Boundary; 2],
    [i, _]: [i32; 2],
    dx: f64,
    _er: f64,
    [_ot, _t]: [f64; 2],
    [_dt, _cdt]: [f64; 2],
    _opt: &(),
    tocomp: ToCompute,
) -> [f64; 1] {
    match tocomp {
        ToCompute::All | ToCompute::Integrated => {
            let ivdx = 1.0 / dx;
            let f = 0;
            let res = ivdx.powf(2.0)
                * (-v[0][boundary[0](i - 2, V)][f] + 16.0 * v[0][boundary[0](i - 1, V)][f]
                    - 30.0 * v[0][boundary[0](i, V)][f]
                    + 16.0 * v[0][boundary[0](i + 1, V)][f]
                    - v[0][boundary[0](i + 2, V)][f])
                / 12.0;
            [res]
        }
        ToCompute::NonIntegrated => [0.0],
    }
}

fn constraints(u: [f64; 1]) -> [f64; 1] {
    u
}

pub fn diffusion() {
    const V: usize = 100;
    let mut vs = [[[0.0]; V]];
    let names = ["u"];
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
        r: ([[1.0]], None),
        dt: 10.0,
        dx: 0.1,
        maxdt: 0.1,
        er: 1e-2,
        t: 0.0,
        tend: 100.0,
        opt: (),
    };
    let (vs, _t, cost, tsteps) = run(context, Integration::FixPointOnly, &names, &constraints);
    println!("cost: {}, tsteps: {}", cost, tsteps);
    let mean_res = vs[0].iter().fold(0.0, |acc, i| acc + i[0]) / V as f64;
    println!("expected mean: {}", mean);
    println!("result mean:   {}", mean_res);
}
