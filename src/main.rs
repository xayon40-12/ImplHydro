use implicit_newton::hydro2d::{f00, f01, f02, hydro2d, Coordinate};

fn gubser(x: f64, y: f64, t: f64) -> [f64; 4] {
    let r = (x * x + y * y).sqrt();
    let e = 2.0f64.powf(8.0 / 3.0)
        / (t.powf(4.0 / 3.0)
            * (1.0 + 2.0 * (t * t + r * r) + (t * t - r * r).powf(2.0)).powf(4.0 / 3.0));
    let k = (2.0 * t * r / (1.0 + t * t + r * r)).atanh();
    let ut = k.cosh();
    let ux;
    let uy;
    if r == 0.0 {
        ux = 0.0;
        uy = 0.0;
    } else {
        ux = x / r * k.sinh();
        uy = y / r * k.sinh();
    }
    [e, ut, ux, uy]
}
fn gubser_err<const V: usize>(v: [[[f64; 4]; V]; V], t: f64, dx: f64) -> [f64; 2] {
    let v2 = ((V - 1) as f64) / 2.0;
    let mut maxerr = 0.0f64;
    let mut meanerr = 0.0f64;
    for i in 0..V {
        for j in 0..V {
            let x = (i as f64 - v2) * dx;
            let y = (j as f64 - v2) * dx;
            let ve = v[j][i][3];
            let [e, _, _, _] = gubser(x, y, t);
            let err = (ve - e).abs() / ve.abs().max(e.abs());
            maxerr = maxerr.max(err);
            meanerr += err;
        }
    }
    meanerr /= (V * V) as f64;
    [maxerr, meanerr]
}

fn main() {
    // diffusion();
    // hydro1d();
    let t0 = 0.6;
    let dx = 0.1;
    // let r = [[1.0]];
    let r = [[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]];
    let (v, t) = hydro2d::<101, 2>(
        0.05,
        1e-3,
        t0,
        t0 + 4.0,
        dx,
        r,
        Coordinate::Milne,
        |x, y| {
            // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
            let [e, ut, ux, uy] = gubser(x, y, t0);
            let vars = [0.0, 0.0, 0.0, e, ut, ux, uy];
            [f00(vars), f01(vars), f02(vars), e]
        },
    );
    let [maxerr, meanerr] = gubser_err(v, t, dx);
    println!("gubser:\nmaxerr: {}\nminerr: {}\n", maxerr, meanerr);
}
