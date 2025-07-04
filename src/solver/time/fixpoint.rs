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
) -> Option<(f64, Box<[[[usize; VX]; VY]; VZ]>, usize)> {
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
    let mut iserr = true;
    let _ko = k.clone();
    let mut fu = k.clone();
    let mut vdtk: Box<[[[[f64; F]; VX]; VY]; VZ]> = boxarray(0.0);
    let mut trdtk: Box<[[[[f64; C]; VX]; VY]; VZ]> = boxarray(0.0);
    // let mut errs: Box<[[[bool; VX]; VY]; VZ]> = boxarray(true);
    let mut nbiter: Box<[[[usize; VX]; VY]; VZ]> = boxarray(0);
    let dt = *dto;
    let mut failed = 1usize;
    let maxfailed = 12;
    let mut maxe = 1e50f64;
    let mut omaxe = maxe;
    let mut alpha = 1.0;
    let mut maxe_count = 0;
    while iserr {
        // let main_time = Instant::now();
        let mut new_coords: Vec<Coord> = Vec::with_capacity(VX * VY * VZ);
        if maxe.is_nan() || maxe > 1e50 || maxe_count > 4 {
            eprintln!("reset {}\t: {}", failed, t);
            failed += 1;
            alpha *= 0.5;
            maxe = 1e50;
            omaxe = maxe;
            // *k = ko.clone();
            const TEST: bool = false;
            if TEST {
                for &Coord { x, y, z } in &coords {
                    for s in 0..S {
                        for f in 0..F {
                            k[s][z][y][x][f] = 0.0;
                            fu[s][z][y][x][f] = 0.0;
                            vdtk[z][y][x][f] = 0.0;
                            trdtk[z][y][x][f] = 0.0;
                        }
                    }
                }
                // coords = gen_coords::<VX, VY, VZ>();
            } else {
                coords = gen_coords::<VX, VY, VZ>();
                *k = boxarray(0.0); // TODO: consider using 'ko' or reset at 0.0 depending on the number of fails
                fu = k.clone();
                *vdtk = [[[[0.0f64; F]; VX]; VY]; VZ];
                *trdtk = [[[[0.0f64; C]; VX]; VY]; VZ];
                // *errs = [[[true; VX]; VY]; VZ];
            }
            if failed >= maxfailed {
                return None;
            }
        }
        let er = err_thr(*t, &vs, &trs);
        let mut maxe = 0.0f64;
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let ct = *t + cdt;
            cfor3d2(
                &coords,
                &mut vdtk,
                &mut trdtk,
                |&Coord { x, y, z }, vdtk, trdtk| {
                    for f in 0..F {
                        vdtk[f] = vs[z][y][x][f];
                        for s1 in 0..S {
                            vdtk[f] += dt * a[s][s1] * k[s1][z][y][x][f];
                        }
                        (*vdtk, *trdtk) = constraints(ct, *vdtk);
                    }
                },
            );
            let (cdt, ovs, otrs) = if cdt == 0.0 {
                (dt, &ovs, &otrs)
            } else {
                (cdt, &vs, &trs)
            };
            cfor3d(&coords, &mut fu[s], |&Coord { x, y, z }, fu| {
                *fu = fun(
                    &k[s],
                    [&ovs, &vdtk],
                    [&otrs, &trdtk],
                    constraints,
                    boundary,
                    [x as i32, y as i32, z as i32],
                    *dxs,
                    [*t, ct],
                    [dt, cdt],
                    opt,
                );
            });
        }
        coords.iter().for_each(|&Coord { z, y, x }| {
            let mut err = false;
            for s in 0..S {
                nbiter[z][y][x] += 1;
                for f in 0..F {
                    let fuf = fu[s][z][y][x][f];
                    let kf = k[s][z][y][x][f];
                    let d = fuf - kf;
                    let e = d.abs();
                    maxe = maxe.max(e);
                    err = err || e.is_nan() || e > er;
                    k[s][z][y][x][f] = kf + alpha * d;
                }
            }
            if err {
                new_coords.push(Coord { x, y, z });
            }
        });
        // let main_elapsed = main_time.elapsed().as_secs_f32();

        // let error_time = Instant::now();
        new_coords.sort_unstable();
        new_coords.dedup();
        if ERROR_PROPAGATION {
            for i in 0..new_coords.len() {
                let Coord { x, y, z } = new_coords[i];
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
                        });
                    }
                }
            }
            new_coords.sort_unstable();
            new_coords.dedup();
        }

        iserr = maxe > er || maxe.is_nan();

        if maxe >= omaxe {
            maxe_count += 1;
        } else {
            maxe_count = 0;
            omaxe = maxe;
        }
        coords = new_coords;
        // let error_elpased = error_time.elapsed().as_secs_f32();
        // let ratio = error_elpased / main_elapsed;
        // if ratio > 5e-2 {
        //     println!(
        //         "error/main: {:e}, error: {:.6}ms, main: {:.6}ms, len: {}",
        //         ratio,
        //         error_elpased,
        //         main_elapsed,
        //         coords.len()
        //     );
        // }
    }
    *ovs = vs.clone();
    *otrs = trs.clone();
    let b = if let Some(b) = b { *b } else { a[S - 1] };
    for s in 0..S {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    for f in 0..F {
                        vs[vz][vy][vx][f] += dt * b[s] * k[s][vz][vy][vx][f];
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
    Some((cost, nbiter, failed - 1))
}
