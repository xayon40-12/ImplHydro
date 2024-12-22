use crate::solver::{
    context::{Arr, Context},
    utils::{cfor3d, cfor3d2, gen_coords, Coord},
    ERROR_PROPAGATION,
};
use boxarray::boxarray;

use super::schemes::Scheme;

pub type ErrThr<
    'a,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> = &'a dyn Fn(f64, &Arr<F, VX, VY, VZ>, &Arr<C, VX, VY, VZ>) -> f64;

pub fn fixpoint<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const S: usize,
>(
    Context {
        fun,
        constraints,
        post_constraints,
        boundary,
        local_interaction,
        vstrs: (vs, trs),
        ovstrs: (ovs, otrs),
        total_diff_vs,
        k,
        r: Scheme { aij: a, bj: b, .. },
        dt: dto,
        dxs,
        maxdt,
        t,
        ot,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
        freezeout_energy: _,
    }: &mut Context<Opt, F, C, VX, VY, VZ, S>,
    err_thr: ErrThr<F, C, VX, VY, VZ>,
) -> Option<(f64, Box<[[[usize; VX]; VY]; VZ]>)> {
    let max_error_increases = 10;

    let b = if let Some(b) = b { *b } else { a[S - 1] };
    let mut coords = gen_coords::<VX, VY, VZ>();
    let [sizex, sizey, sizez] = *local_interaction;
    let dxyz: Vec<[i32; 3]> = {
        let mut d = vec![];
        for dx in -sizex..=sizex {
            d.push([dx, 0, 0]);
        }
        for dy in 1..=sizey {
            d.push([0, -dy, 0]);
            d.push([0, dy, 0]);
        }
        for dz in 1..=sizez {
            d.push([0, 0, -dz]);
            d.push([0, 0, dz]);
        }
        d
    };
    *dto = maxdt.min(*dto);
    let mut ko = k.clone();
    let mut fu = k.clone();
    let mut vdtk: Box<[[[[f64; F]; VX]; VY]; VZ]> = boxarray(0.0);
    let mut trdtk: Box<[[[[f64; C]; VX]; VY]; VZ]> = boxarray(0.0);
    // let mut errs: Box<[[[bool; VX]; VY]; VZ]> = boxarray(true);
    let mut nbiter: Box<[[[usize; VX]; VY]; VZ]> = boxarray(0);
    let mut dts: Box<[[[f64; VX]; VY]; VZ]> = boxarray(*dto);
    let mut ts: Box<[[[f64; VX]; VY]; VZ]> = boxarray(*t);
    while coords.len() > 0 {
        let mut new_coords: Vec<Coord> = Vec::with_capacity(VX * VY * VZ);
        let er = err_thr(*t, &vs, &trs);
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            cfor3d2(
                &coords,
                &mut vdtk,
                &mut trdtk,
                |&Coord { x, y, z, .. }, vdtk, trdtk| {
                    let t = ts[z][y][x];
                    let dt = dts[z][y][x];
                    let ct = t + c * dt;
                    for f in 0..F {
                        vdtk[f] = vs[z][y][x][f];
                        for s1 in 0..S {
                            vdtk[f] += dt * a[s][s1] * k[s1][z][y][x][f];
                        }
                        (*vdtk, *trdtk) = constraints(ct, *vdtk);
                    }
                },
            );
            let (ovs, otrs) = if c == 0.0 { (&ovs, &otrs) } else { (&vs, &trs) };
            cfor3d(&coords, &mut fu[s], |&Coord { x, y, z, .. }, fu| {
                let t = ts[z][y][x];
                let dt = dts[z][y][x];
                let cdt = if c == 0.0 { dt } else { c * dt };
                let ct = t + c * dt;
                *fu = fun(
                    &k[s],
                    [&ovs, &vdtk],
                    [&otrs, &trdtk],
                    constraints,
                    boundary,
                    [x as i32, y as i32, z as i32],
                    *dxs,
                    [t, ct],
                    [dt, cdt],
                    opt,
                );
            });
        }
        coords.iter().for_each(
            |&Coord {
                 z,
                 y,
                 x,
                 mut remaining,
                 mut error_increases,
                 mut max_err,
             }| {
                let mut is_err = false;
                let mut current_max_err = 0.0f64;
                for s in 0..S {
                    nbiter[z][y][x] += 1;
                    for f in 0..F {
                        let fuf = fu[s][z][y][x][f];
                        let kf = k[s][z][y][x][f];
                        let d = fuf - kf;
                        let e = d.abs();
                        current_max_err = current_max_err.max(e);
                        is_err |= e.is_nan() || e > er; // * dts[z][y][x] / *dto;
                        k[s][z][y][x][f] = kf + d;
                    }
                }
                let error_increased = current_max_err >= max_err;
                max_err = current_max_err;
                if is_err {
                    if error_increased {
                        error_increases += 1;
                        if error_increases > max_error_increases {
                            remaining *= 2;
                            dts[z][y][x] *= 0.5;
                            for s in 0..S {
                                for f in 0..F {
                                    k[s][z][y][x][f] = ko[s][z][y][x][f];
                                }
                            }
                            max_err = f64::MAX;
                        }
                    } else {
                        error_increases = 0;
                    }
                    new_coords.push(Coord {
                        x,
                        y,
                        z,
                        remaining,
                        error_increases,
                        max_err,
                    });
                } else {
                    if remaining > 1 {
                        remaining -= 1;
                        for s in 0..S {
                            for f in 0..F {
                                vs[z][y][x][f] += dts[z][y][x] * b[s] * k[s][z][y][x][f];
                                ko[s][z][y][x][f] = k[s][z][y][x][f];
                            }
                        }
                        ts[z][y][x] += dts[z][y][x];
                        let t = ts[z][y][x];
                        (vs[z][y][x], trs[z][y][x]) = constraints(t, vs[z][y][x]);
                        if let Some(post) = post_constraints {
                            (vs[z][y][x], trs[z][y][x]) = post(t, vs[z][y][x]);
                        }
                        new_coords.push(Coord {
                            x,
                            y,
                            z,
                            remaining,
                            error_increases: 0,
                            max_err,
                        });
                    }
                }
            },
        );

        new_coords.sort_by_key(|&Coord { z, y, x, .. }| (z, y, x));
        const DEBUG: bool = false;
        if DEBUG {
            println!("coords.len(): {}", coords.len());
            if coords.len() == 1 {
                let Coord {
                    z,
                    y,
                    x,
                    remaining,
                    error_increases,
                    max_err,
                } = coords[0];
                println!("k: {:?}", k[0][z][y][x]);
                println!("dt: {:e}", dts[z][y][x]);
                println!("remaining: {remaining}");
                println!("max_err: {max_err:e}");
                if dts[z][y][x] < 1e-10 || k[0][z][y][x][0].is_nan() {
                    panic!("dt null or nan");
                }
            }
        }
        if ERROR_PROPAGATION {
            for i in 0..new_coords.len() {
                let c = new_coords[i];
                let Coord { x, y, z, .. } = c;
                for &[dx, dy, dz] in &dxyz {
                    let i = x as i32 + dx;
                    let j = y as i32 + dy;
                    let k = z as i32 + dz;
                    if i >= 0 && i < VX as i32 && j >= 0 && j < VY as i32 && k >= 0 && k < VZ as i32
                    {
                        new_coords.push(Coord {
                            x: i as usize,
                            y: j as usize,
                            z: k as usize,
                            ..c
                        });
                    }
                }
            }
            new_coords.sort_by_key(|&Coord { z, y, x, .. }| (z, y, x));
            new_coords.dedup_by_key(|&mut Coord { z, y, x, .. }| (z, y, x));
        }

        coords = new_coords;
    }
    *ovs = vs.clone();
    *otrs = trs.clone();
    for s in 0..S {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    for f in 0..F {
                        vs[vz][vy][vx][f] += dts[vz][vy][vx] * b[s] * k[s][vz][vy][vx][f];
                    }
                }
            }
        }
    }
    let cost = (nbiter
        .iter()
        .flat_map(|e| e.iter())
        .flat_map(|e| e.iter())
        .map(|v| *v)
        .reduce(|acc, a| acc + a)
        .unwrap()) as f64
        / (VX * VY * VZ) as f64;
    *ot = *t;
    let dt = *dto;
    *t += dt;
    *dto = maxdt.min(dt * 1.1);
    for vz in 0..VZ {
        for vy in 0..VY {
            for vx in 0..VX {
                let tmp = vs[vz][vy][vx];
                (vs[vz][vy][vx], trs[vz][vy][vx]) = constraints(*t, vs[vz][vy][vx]);
                for f in 0..F {
                    total_diff_vs[vz][vy][vx][f] += (tmp[f] - vs[vz][vy][vx][f]).abs();
                }
            }
        }
    }
    if let Some(post) = post_constraints {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    let tmp = vs[vz][vy][vx];
                    (vs[vz][vy][vx], trs[vz][vy][vx]) = post(*t, vs[vz][vy][vx]);
                    for f in 0..F {
                        total_diff_vs[vz][vy][vx][f] += (tmp[f] - vs[vz][vy][vx][f]).abs();
                    }
                }
            }
        }
    }
    Some((cost, nbiter))
}
