use rayon::prelude::*;

pub fn zero(_j: i32, _n: usize) -> usize {
    0
}

pub fn zeros<const VY: usize, const VX: usize, const F: usize>(
    _: &[[[f64; F]; VX]; VY],
) -> [[[f64; F]; VX]; VY] {
    [[[0.0; F]; VX]; VY]
}

pub fn _periodic(j: i32, n: usize) -> usize {
    if j >= n as i32 {
        j as usize - n
    } else if j < 0 {
        (n as i32 + j) as usize
    } else {
        j as usize
    }
}
pub fn periodic<const F: usize, const C: usize, const VX: usize, const VY: usize>(
    [x, y]: [i32; 2],
    vs: [[[f64; F]; VX]; VY],
    trs: [[[f64; C]; VX]; VY],
) -> ([f64; F], [f64; C]) {
    let j = _periodic(y, VY);
    let i = _periodic(x, VX);
    (vs[j][i], trs[j][i])
}

pub fn _ghost(j: i32, n: usize) -> usize {
    if j >= n as i32 {
        n - 1
    } else if j < 0 {
        0
    } else {
        j as usize
    }
}
pub fn ghost<const F: usize, const VX: usize, const VY: usize>(
    [x, y]: [i32; 2],
    vs: &[[[f64; F]; VX]; VY],
) -> [f64; F] {
    let j = _ghost(y, VY);
    let i = _ghost(x, VX);
    vs[j][i]
}

pub fn _cubical<const F: usize>(x: i32, vs: [[f64; F]; 4]) -> [f64; F] {
    let x = x as f64;
    let mut res = [0.0f64; F];
    let mat = [
        [-1.0 / 6.0, 1.0 / 2.0, -1.0 / 2.0, 1.0 / 6.0],
        [1.0, -5.0 / 2.0, 2.0, -1.0 / 2.0],
        [-11.0 / 6.0, 3.0, -3.0 / 2.0, 1.0 / 3.0],
        [1.0, 0.0, 0.0, 0.0],
    ];
    for f in 0..F {
        let mut m = [0.0; 4];
        for j in 0..4 {
            for i in 0..4 {
                m[j] += mat[j][i] * vs[i][f];
            }
        }

        res[f] = m[0] * x.powi(3) + m[1] * x.powi(2) + m[2] * x + m[3];
    }

    res
}

pub fn cubical<const F: usize, const VX: usize, const VY: usize>(
    [x, y]: [i32; 2],
    vs: &[[[f64; F]; VX]; VY],
) -> [f64; F] {
    if x < 0 || x >= VX as i32 {
        let (s, p) = if x < 0 {
            (0, x)
        } else {
            (VX - 5, x - VX as i32 + 4)
        };
        let mut cvs = [[0.0f64; F]; 4];
        let j = y as usize;
        for i in 0..4 {
            cvs[i] = vs[j][s + i];
        }
        _cubical(p, cvs)
    } else if y < 0 || y >= VY as i32 {
        let (s, p) = if y < 0 {
            (0, y)
        } else {
            (VY - 5, y - VY as i32 + 4)
        };
        let mut cvs = [[0.0f64; F]; 4];
        let i = x as usize;
        for j in 0..4 {
            cvs[j] = vs[s + j][i];
        }
        _cubical(p, cvs)
    } else {
        let j = y as usize;
        let i = x as usize;
        vs[j][i]
    }
}

pub fn flux_limiter(theta: f64, a: f64, b: f64, c: f64) -> f64 {
    let minmod2 = |a: f64, b: f64| (a.signum() + b.signum()) / 2.0 * a.abs().min(b.abs());
    minmod2(theta * (b - a), minmod2((c - a) / 2.0, theta * (c - b)))
}

pub struct Coord {
    pub x: usize,
    pub y: usize,
}

pub fn pfor2d<T: Send, const VX: usize, const VY: usize>(
    vss: &mut [[T; VX]; VY],
    f: &(dyn Fn((Coord, &mut T)) + Sync),
) {
    vss.par_iter_mut()
        .enumerate()
        .flat_map(|(vy, vs)| {
            vs.par_iter_mut()
                .enumerate()
                .map(move |(vx, v)| (Coord { x: vx, y: vy }, v))
        })
        .for_each(f);
}
pub fn pfor2d2<T: Send, U: Send, const VX: usize, const VY: usize>(
    vss: &mut [[T; VX]; VY],
    wss: &mut [[U; VX]; VY],
    f: &(dyn Fn((Coord, &mut T, &mut U)) + Sync),
) {
    vss.par_iter_mut()
        .zip(wss.par_iter_mut())
        .enumerate()
        .flat_map(|(vy, (vs, ws))| {
            vs.par_iter_mut()
                .zip(ws.par_iter_mut())
                .enumerate()
                .map(move |(vx, (v, w))| (Coord { x: vx, y: vy }, v, w))
        })
        .for_each(f);
}
