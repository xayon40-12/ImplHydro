use crate::FLOAT;

pub fn newton(
    mut v: FLOAT,
    f: impl Fn(FLOAT) -> FLOAT,
    constraint: impl Fn(FLOAT) -> FLOAT,
) -> FLOAT {
    let er = (10.0 as FLOAT).powi(-(FLOAT::DIGITS as i32) * 2 / 3);
    let mut fv = f(v);
    let mut e = fv.abs();
    let mut i = 0;
    let maxi = 1000;
    let mut ov;
    while e >= er && i < maxi {
        i += 1;
        let ff = (f(v + er) - fv) / er;
        ov = v;
        v -= fv / (er + ff);
        v = constraint(v);
        fv = f(v);
        e = fv.abs() + (ov - v).abs() / (ov.abs().max(v.abs()) + er);
    }
    if i == maxi {
        FLOAT::NAN
    } else {
        v
    }
}
