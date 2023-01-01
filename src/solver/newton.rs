pub fn newton(er: f64, mut v: f64, f: impl Fn(f64) -> f64) -> f64 {
    let mut e = er + 1.0;
    while e >= er {
        let fv = f(v);
        let ff = (f(v + er) - fv) / er;
        v -= fv / ff;
        e = fv.abs();
    }
    v
}
