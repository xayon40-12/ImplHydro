use boxarray::boxarray;

use super::context::{Arr, BArr, DIM};

pub fn zero(_j: i32, _n: usize) -> usize {
    0
}

pub fn zeros<const F: usize, const VX: usize, const VY: usize, const VZ: usize>(
) -> BArr<F, VX, VY, VZ> {
    boxarray(0.0)
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
pub fn periodic<const F: usize, const VX: usize, const VY: usize, const VZ: usize>(
    [x, y, z]: [i32; DIM],
    vs: &Arr<F, VX, VY, VZ>,
) -> [f64; F] {
    let k = _periodic(z, VZ);
    let j = _periodic(y, VY);
    let i = _periodic(x, VX);
    vs[k][j][i]
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
pub fn ghost<const F: usize, const VX: usize, const VY: usize, const VZ: usize>(
    [x, y, z]: [i32; DIM],
    vs: &Arr<F, VX, VY, VZ>,
) -> [f64; F] {
    let k = _ghost(z, VZ);
    let j = _ghost(y, VY);
    let i = _ghost(x, VX);
    vs[k][j][i]
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

pub fn cubical<const F: usize, const VX: usize, const VY: usize, const VZ: usize>(
    [x, y, z]: [i32; DIM],
    vs: &Arr<F, VX, VY, VZ>,
) -> [f64; F] {
    let k = z as usize;
    let j = y as usize;
    let i = x as usize;
    if x < 0 || x >= VX as i32 {
        let (s, p) = if x < 0 {
            (0, x)
        } else {
            (VX - 5, x - VX as i32 + 4)
        };
        let mut cvs = [[0.0f64; F]; 4];
        for i in 0..4 {
            cvs[i] = vs[k][j][s + i];
        }
        _cubical(p, cvs)
    } else if y < 0 || y >= VY as i32 {
        let (s, p) = if y < 0 {
            (0, y)
        } else {
            (VY - 5, y - VY as i32 + 4)
        };
        let mut cvs = [[0.0f64; F]; 4];
        for j in 0..4 {
            cvs[j] = vs[k][s + j][i];
        }
        _cubical(p, cvs)
    } else if z < 0 || z >= VZ as i32 {
        let (s, p) = if z < 0 {
            (0, z)
        } else {
            (VZ - 5, z - VZ as i32 + 4)
        };
        let mut cvs = [[0.0f64; F]; 4];
        for k in 0..4 {
            cvs[k] = vs[s + k][j][i];
        }
        _cubical(p, cvs)
    } else {
        vs[k][j][i]
    }
}

pub fn flux_limiter(theta: f64, a: f64, b: f64, c: f64) -> f64 {
    let minmod2 = |a: f64, b: f64| (a.signum() + b.signum()) / 2.0 * a.abs().min(b.abs());
    minmod2(theta * (b - a), minmod2((c - a) / 2.0, theta * (c - b)))
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
#[repr(C)]
pub struct Coord {
    pub z: usize,
    pub y: usize,
    pub x: usize,
    pub remaining: usize,
    pub error_increases: usize,
    pub max_err: f64,
}

pub fn gen_coords<const VX: usize, const VY: usize, const VZ: usize>() -> Vec<Coord> {
    let mut coords = Vec::with_capacity(VX * VY * VZ);
    for z in 0..VZ {
        for y in 0..VY {
            for x in 0..VX {
                coords.push(Coord {
                    x,
                    y,
                    z,
                    remaining: 1,
                    error_increases: 0,
                    max_err: f64::MAX,
                });
            }
        }
    }
    coords
}

macro_rules! gen_cfor3d {
    ($n: ident, $($t: ident: $T: ident),*) => {
        pub fn $n<$($T,)* const VX: usize, const VY: usize, const VZ: usize>(
            coords: &[Coord],
            $($t: &mut [[[$T; VX]; VY]; VZ],)*
            mut f: impl FnMut(&Coord, $(&mut $T),*),
        ) {
            coords.iter().for_each(|c| f(c, $(&mut $t[c.z][c.y][c.x]),*));
        }
    };
}
gen_cfor3d!(cfor3d, t1: T1);
gen_cfor3d!(cfor3d2, t1: T1, t2: T2);
gen_cfor3d!(cfor3d3, t1: T1, t2: T2, t3: T3);
gen_cfor3d!(cfor3d4, t1: T1, t2: T2, t3: T3, t4: T4);
gen_cfor3d!(cfor3d5, t1: T1, t2: T2, t3: T3, t4: T4, t5: T5);
