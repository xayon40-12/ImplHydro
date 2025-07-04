use std::{mem::MaybeUninit, ptr::addr_of_mut};

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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
#[repr(C)]
pub struct Coord {
    pub z: usize,
    pub y: usize,
    pub x: usize,
}

pub fn gen_coords<const VX: usize, const VY: usize, const VZ: usize>() -> Vec<Coord> {
    let mut coords = Vec::with_capacity(VX * VY * VZ);
    for z in 0..VZ {
        for y in 0..VY {
            for x in 0..VX {
                coords.push(Coord { x, y, z });
            }
        }
    }
    coords
}

pub fn cfor3dn<T: Send, const N: usize, const VX: usize, const VY: usize, const VZ: usize>(
    coords: &[Coord],
    mut nsss: [&mut [[[T; VX]; VY]; VZ]; N],
    f: impl Fn(&Coord, [&mut T; N]) + Sync,
) {
    coords.iter().for_each(|c| {
        let tsss: [&mut T; N] = unsafe {
            let mut uninit: MaybeUninit<[&mut T; N]> = MaybeUninit::uninit();
            let ptr = uninit.as_mut_ptr();
            nsss.iter_mut().enumerate().for_each(|(i, sss)| {
                addr_of_mut!((*ptr)[i]).write(&mut sss[c.z][c.y][c.x]);
            });
            uninit.assume_init()
        };
        f(c, tsss)
    });
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

pub fn pfor3d<T: Send, const VX: usize, const VY: usize, const VZ: usize>(
    vsss: &mut [[[T; VX]; VY]; VZ],
    f: impl Fn((Coord, &mut T)) + Sync,
) {
    vsss.iter_mut()
        .enumerate()
        .flat_map(|(vz, vss)| {
            vss.iter_mut().enumerate().flat_map(move |(vy, vs)| {
                vs.iter_mut().enumerate().map(move |(vx, v)| {
                    (
                        Coord {
                            x: vx,
                            y: vy,
                            z: vz,
                        },
                        v,
                    )
                })
            })
        })
        .for_each(f);
}
pub fn pfor3d2<T: Send, U: Send, const VX: usize, const VY: usize, const VZ: usize>(
    vsss: &mut [[[T; VX]; VY]; VZ],
    wsss: &mut [[[U; VX]; VY]; VZ],
    f: impl Fn((Coord, &mut T, &mut U)) + Sync,
) {
    vsss.iter_mut()
        .zip(wsss.iter_mut())
        .enumerate()
        .flat_map(|(vz, (vss, wss))| {
            vss.iter_mut()
                .zip(wss.iter_mut())
                .enumerate()
                .flat_map(move |(vy, (vs, ws))| {
                    vs.iter_mut()
                        .zip(ws.iter_mut())
                        .enumerate()
                        .map(move |(vx, (v, w))| {
                            (
                                Coord {
                                    x: vx,
                                    y: vy,
                                    z: vz,
                                },
                                v,
                                w,
                            )
                        })
                })
        })
        .for_each(f);
}
pub fn pfor3d3<U: Send, V: Send, W: Send, const VX: usize, const VY: usize, const VZ: usize>(
    usss: &mut [[[U; VX]; VY]; VZ],
    vsss: &mut [[[V; VX]; VY]; VZ],
    wsss: &mut [[[W; VX]; VY]; VZ],
    f: impl Fn((Coord, &mut U, &mut V, &mut W)) + Sync,
) {
    usss.iter_mut()
        .zip(vsss.iter_mut())
        .zip(wsss.iter_mut())
        .enumerate()
        .flat_map(|(vz, ((uss, vss), wss))| {
            uss.iter_mut()
                .zip(vss.iter_mut())
                .zip(wss.iter_mut())
                .enumerate()
                .flat_map(move |(vy, ((us, vs), ws))| {
                    us.iter_mut()
                        .zip(vs.iter_mut())
                        .zip(ws.iter_mut())
                        .enumerate()
                        .map(move |(vx, ((u, v), w))| {
                            (
                                Coord {
                                    x: vx,
                                    y: vy,
                                    z: vz,
                                },
                                u,
                                v,
                                w,
                            )
                        })
                })
        })
        .for_each(f);
}
pub fn pfor3d4<
    T1: Send,
    T2: Send,
    T3: Send,
    T4: Send,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
>(
    t1sss: &mut [[[T1; VX]; VY]; VZ],
    t2sss: &mut [[[T2; VX]; VY]; VZ],
    t3sss: &mut [[[T3; VX]; VY]; VZ],
    t4sss: &mut [[[T4; VX]; VY]; VZ],
    f: impl Fn((Coord, &mut T1, &mut T2, &mut T3, &mut T4)) + Sync,
) {
    t1sss
        .iter_mut()
        .zip(t2sss.iter_mut())
        .zip(t3sss.iter_mut())
        .zip(t4sss.iter_mut())
        .enumerate()
        .flat_map(|(vz, (((t1ss, t2ss), t3ss), t4ss))| {
            t1ss.iter_mut()
                .zip(t2ss.iter_mut())
                .zip(t3ss.iter_mut())
                .zip(t4ss.iter_mut())
                .enumerate()
                .flat_map(move |(vy, (((t1s, t2s), t3s), t4s))| {
                    t1s.iter_mut()
                        .zip(t2s.iter_mut())
                        .zip(t3s.iter_mut())
                        .zip(t4s.iter_mut())
                        .enumerate()
                        .map(move |(vx, (((t1, t2), t3), t4))| {
                            (
                                Coord {
                                    x: vx,
                                    y: vy,
                                    z: vz,
                                },
                                t1,
                                t2,
                                t3,
                                t4,
                            )
                        })
                })
        })
        .for_each(f);
}
