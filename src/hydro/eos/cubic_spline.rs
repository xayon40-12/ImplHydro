use crate::boxarray;

pub fn cubic(x: f64, p: &[f64; 4], _a: f64) -> f64 {
    p[0] + x * (p[1] + x * (p[2] + x * p[3]))
}

pub fn cubic_diff(x: f64, p: &[f64; 4], a: f64) -> f64 {
    (p[1] + x * (2.0 * p[2] + x * 3.0 * p[3])) / a
}

pub fn cubic_diff2(x: f64, p: &[f64; 4], a: f64) -> f64 {
    (2.0 * p[2] + x * 6.0 * p[3]) / (a * a)
}

pub fn cubic_spline_0<const N: usize, const M: usize>(
    idx: usize,
    id: usize,
    arr: &[[f64; M]; N],
) -> Box<[[f64; 4]; N]> {
    let x0 = arr[0][idx];
    let x1 = arr[1][idx];
    let x1x0 = x1 - x0;
    let a = x0 / x1x0;
    let aa = a * a;

    let mut tmp = vec![[0.0f64; 3]; N + 1];
    // WARNING: as the first segment down to 0 is not the same lenght as the other segments, the spline does not work properly for this first segment
    tmp[0][0] = 0.0;
    tmp[0][2] = 2.0 * a;
    for i in 0..N {
        tmp[i + 1][0] = arr[i][id];
        tmp[i + 1][2] = 4.0;
    }
    let m01 = a;
    let m10 = a;
    tmp[1][2] = 2.0 * a * (1.0 + a);
    let m12 = aa;
    tmp[N][2] = 2.0;

    tmp[0][1] = 3.0 * (tmp[1][0] - tmp[0][0]);
    tmp[1][1] = 3.0 * (aa * tmp[2][0] + (1.0 - aa) * tmp[1][0] - tmp[0][0]);
    tmp[N][1] = 3.0 * (tmp[N][0] - tmp[N - 1][0]);
    for i in 2..N {
        tmp[i][1] = 3.0 * (tmp[i + 1][0] - tmp[i - 1][0]);
    }

    let d = m10 / tmp[0][2];
    tmp[1][2] -= d * m01;
    tmp[1][1] -= d * tmp[0][1];
    let d = 1.0 / tmp[1][2];
    tmp[2][2] -= d * m12;
    tmp[2][1] -= d * tmp[1][1];
    for i in 3..=N {
        let d = 1.0 / tmp[i - 1][2];
        tmp[i][2] -= d;
        tmp[i][1] -= d * tmp[i - 1][1];
    }
    for i in (2..N).rev() {
        tmp[i + 1][1] /= tmp[i + 1][2];
        tmp[i][1] -= tmp[i + 1][1];
    }
    tmp[2][1] /= tmp[2][2];
    tmp[1][1] -= tmp[2][1] * m12;
    tmp[1][1] /= tmp[1][2];
    tmp[0][1] -= tmp[1][1] * m01;
    tmp[0][1] /= tmp[0][2];

    let mut pols: Box<[[f64; 4]; N]> = boxarray(0.0f64);
    pols[0][0] = tmp[0][0];
    pols[0][1] = a * tmp[0][1];
    pols[0][2] = 3.0 * (tmp[1][0] - tmp[0][0]) - a * (2.0 * tmp[0][1] + tmp[1][1]);
    pols[0][3] = 2.0 * (tmp[0][0] - tmp[1][0]) + a * (tmp[0][1] + tmp[1][1]);
    for i in 1..N {
        pols[i][0] = tmp[i][0];
        pols[i][1] = tmp[i][1];
        pols[i][2] = 3.0 * (tmp[i + 1][0] - tmp[i][0]) - 2.0 * tmp[i][1] - tmp[i + 1][1];
        pols[i][3] = 2.0 * (tmp[i][0] - tmp[i + 1][0]) + tmp[i][1] + tmp[i + 1][1];
    }

    pols
}
