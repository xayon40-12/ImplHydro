use rayon::prelude::*;

pub fn zero(_j: i32, _n: usize) -> usize {
    0
}

pub fn periodic(j: i32, n: usize) -> usize {
    if j >= n as i32 {
        j as usize - n
    } else if j < 0 {
        (n as i32 + j) as usize
    } else {
        j as usize
    }
}

pub fn ghost(j: i32, n: usize) -> usize {
    if j >= n as i32 {
        n - 1
    } else if j < 0 {
        0
    } else {
        j as usize
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
