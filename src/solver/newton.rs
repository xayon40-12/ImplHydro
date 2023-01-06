pub fn newton(er: f64, mut v: f64, f: impl Fn(f64) -> f64) -> f64 {
    let mut e = er + 1.0;
    let mut i = 0;
    let maxi = 1000;
    while e >= er && i < maxi {
        i += 1;
        let fv = f(v);
        let ff = (f(v + er) - fv) / er;
        v -= fv / ff;
        e = fv.abs();
    }
    if i == maxi {
        f64::NAN
    } else {
        v
    }
}
