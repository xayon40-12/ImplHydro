pub fn newton(er: f64, mut v: f64, f: impl Fn(f64) -> f64, constraint: impl Fn(f64) -> f64) -> f64 {
    let mut fv = f(v);
    let mut e = fv.abs();
    let mut i = 0;
    let maxi = 100;
    while e >= er && i < maxi {
        i += 1;
        let ff = (f(v + er) - fv) / er;
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
