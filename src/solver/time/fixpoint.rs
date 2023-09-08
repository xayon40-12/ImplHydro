use rayon::prelude::*;

use crate::{
    boxarray::boxarray,
    solver::{
        context::{Arr, Context},
        utils::{pfor3d, pfor3d2, Coord},
    },
};

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
        dx,
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
    let [sizex, sizey, sizez] = *local_interaction;
    *dto = maxdt.min(*dto);
    let mut iserr = true;
    let ko = k.clone();
    let mut fu = k.clone();
    let mut vdtk: Box<[[[[f64; F]; VX]; VY]; VZ]> = boxarray(0.0);
    let mut trdtk: Box<[[[[f64; C]; VX]; VY]; VZ]> = boxarray(0.0);
    let mut errs: Box<[[[bool; VX]; VY]; VZ]> = boxarray(true);
    let mut nbiter: Box<[[[usize; VX]; VY]; VZ]> = boxarray(0);
    let dt = *dto;
    let mut failed = 1usize;
    let maxfailed = 10;
    let mut maxe = 1e50f64;
    let mut omaxe = maxe;
    let mut alpha = 1.0;
    let mut maxe_count = 0;
    while iserr {
        if maxe > 1e50 || maxe_count > 4 {
            eprintln!("fail {}", failed);
            failed += 1;
            alpha *= (*dx / dt).min(0.5);
            maxe = 1e50;
            omaxe = maxe;
            *k = ko.clone();
            fu = k.clone();
            *vdtk = [[[[0.0f64; F]; VX]; VY]; VZ];
            *trdtk = [[[[0.0f64; C]; VX]; VY]; VZ];
            *errs = [[[true; VX]; VY]; VZ];
            if failed >= maxfailed {
                return None;
            }
        }
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let ct = *t + cdt;
            pfor3d2(&mut vdtk, &mut trdtk, &|(
                Coord { x, y, z },
                vdtk,
                trdtk,
            )| {
                for f in 0..F {
                    vdtk[f] = vs[z][y][x][f];
                    for s1 in 0..S {
                        vdtk[f] += dt * a[s][s1] * k[s1][z][y][x][f];
                    }
                    (*vdtk, *trdtk) = constraints(ct, *vdtk);
                }
            });
            let (cdt, ovs, otrs) = if cdt == 0.0 {
                (dt, &ovs, &otrs)
            } else {
                (cdt, &vs, &trs)
            };
            pfor3d(&mut fu[s], &|(Coord { x, y, z }, fu)| {
                if errs[z][y][x] {
                    *fu = fun(
                        [&ovs, &vdtk],
                        [&otrs, &trdtk],
                        constraints,
                        boundary,
                        [x as i32, y as i32, z as i32],
                        *dx,
                        [*t, ct],
                        [dt, cdt],
                        opt,
                    );
                }
            });
        }

        let er = err_thr(*t, &vs, &trs);
        pfor3d2(&mut errs, &mut nbiter, &|(
            Coord { x, y, z },
            errs,
            nbiter,
        )| {
            if *errs {
                *errs = false;
                *nbiter += 1;
                for s in 0..S {
                    for f in 0..F {
                        let fuf = fu[s][z][y][x][f];
                        let kf = k[s][z][y][x][f];
                        let e = (fuf - kf).abs();
                        if e.is_nan() {
                            panic!(
                                "\n\nNaN encountered in fixpoint iteration.\nfields: {:?}\n\n",
                                fu[s][z][y][x]
                            );
                        }
                        *errs |= e > er || e.is_nan();
                    }
                }
            }
        });

        iserr = errs
            .par_iter()
            .flat_map(|e| e.par_iter())
            .flat_map(|e| e.par_iter())
            .map(|v| *v)
            .reduce(|| false, |acc, a| acc || a);

        maxe = 0.0;
        let mut errr: Box<[[[f64; VX]; VY]; VZ]> = boxarray(0.0);
        for s in 0..S {
            for z in 0..VZ {
                for y in 0..VY {
                    for x in 0..VX {
                        for f in 0..F {
                            let kk = k[s][z][y][x][f];
                            let ff = fu[s][z][y][x][f];
                            let d = ff - kk;
                            errr[z][y][x] = errr[z][y][x].max(d.abs());
                            maxe = maxe.max(d.abs());
                            k[s][z][y][x][f] = kk + alpha * d;
                        }
                    }
                }
            }
        }
        if maxe >= omaxe {
            maxe_count += 1;
        } else {
            maxe_count = 0;
            omaxe = maxe;
        }

        let debug = false;
        // let debug = true;
        if debug {
            let miter = nbiter
                .iter()
                .flat_map(|v| v.iter())
                .flat_map(|v| v.iter())
                .max()
                .unwrap();
            println!(
                "t: {:.3e}, ffka: {:.3e}, maxe: {:.3e}, omaxe: {:.3e}, lim: {:.3e}, miter: {}",
                t, alpha, maxe, omaxe, er, miter
            );

            let col: Vec<char> = " .:-=+*%#@".chars().collect();
            let mc = col.len() - 1;
            let iter = nbiter[0]
                .iter()
                .map(|v| {
                    v.iter()
                        .map(|v| format!("{s}{s}", s = col[v * mc / miter].to_string()))
                        .collect::<Vec<_>>()
                        .join("")
                })
                .collect::<Vec<_>>();
            let err = errr[0]
                .iter()
                .map(|v| {
                    v.iter()
                        .map(|v| {
                            format!(
                                "{s}{s}",
                                s = col[(((v - er) * (mc - 1) as f64) / (maxe - er).max(1e-100)
                                    + 1.0)
                                    .max(0.0) as usize]
                                    .to_string()
                            )
                        })
                        .collect::<Vec<_>>()
                        .join("")
                })
                .collect::<Vec<_>>();

            println!(
                "{}",
                iter.iter()
                    .zip(err.iter())
                    .map(|(a, b)| format!("{}{}", a, b))
                    .collect::<Vec<_>>()
                    .join("\n")
            );
        }

        let mut tmperrs: Box<[[[bool; VX]; VY]; VZ]> = boxarray(false);
        pfor3d(&mut tmperrs, &|(Coord { x, y, z }, tmperrs)| {
            let mut update = |dx, dy, dz| {
                let i = x as i32 + dx;
                let j = y as i32 + dy;
                let k = z as i32 + dz;
                if i >= 0 && i < VX as i32 && j >= 0 && j < VY as i32 && k >= 0 && k < VZ as i32 {
                    *tmperrs |= errs[k as usize][j as usize][i as usize];
                }
            };
            for dz in -sizez..=sizez {
                if dz == 0 {
                    for dy in -sizey..=sizey {
                        if dy == 0 {
                            for dx in -sizex..=sizex {
                                update(dx, 0, 0);
                            }
                        } else {
                            update(0, dy, 0);
                        }
                    }
                } else {
                    update(0, 0, dz);
                }
            }
        });
        errs = tmperrs;
    }
    *ovs = vs.clone();
    *otrs = trs.clone();
    let b = if let Some(b) = b { *b } else { a[S - 1] };
    for f in 0..F {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    for s in 0..S {
                        vs[vz][vy][vx][f] += dt * b[s] * k[s][vz][vy][vx][f];
                    }
                }
            }
        }
    }
    let cost = (S * nbiter
        .par_iter()
        .flat_map(|e| e.par_iter())
        .flat_map(|e| e.par_iter())
        .map(|v| *v)
        .reduce(|| 0, |acc, a| acc + a)) as f64;
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
    Some((cost / (VX * VY * VZ) as f64, nbiter, failed - 1))
}
