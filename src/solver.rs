use crate::newton::{newton, Boundary, Fun};

pub fn run<'a, const F: usize, const V: usize, const S: usize>(
    f: (Fun<'a, F, V>, Boundary<'a>, i32),
    mut vs: [[f64; V]; F],
    r: [[f64; S]; S],
    dt: f64,
    tbeg: f64,
    tend: f64,
) {
    let mul = f.2 * 2 + 1;
    let tsteps = ((tend - tbeg) / dt) as usize;
    let mut tot_f = 0;
    for i in 0..tsteps {
        let (nvs, c) = newton(f, vs, dt, &r);
        let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        vs = nvs;
        tot_f += c;
        if i % (1 + tsteps / 10) == 0 {
            println!("i: {}, c: {}", i, c);
        }
    }
    println!("v: {:?}", vs);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
}
