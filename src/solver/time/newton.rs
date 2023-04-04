pub fn newton(er: f64, mut v: f64, f: impl Fn(f64) -> f64, constraint: impl Fn(f64) -> f64) -> f64 {
    let mut fv = f(v);
    let mut e = fv.abs();
    let mut i = 0;
    let maxi = 100;
    let mut ov = v;
    // while (e >= er || (ov - v).abs() / (ov.abs().max(v.abs()) + 1e-15) >= er) && i < maxi {
    while e >= er && i < maxi {
        i += 1;
        let ff = (f(v + er) - fv) / er;
        ov = v;
        v -= fv / ff;
        v = constraint(v);
        fv = f(v);
        e = fv.abs();
    }
    if i == maxi {
        f64::NAN
    } else {
        v
    }
}
