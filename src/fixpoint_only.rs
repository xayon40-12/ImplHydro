use rayon::prelude::*;

use crate::context::{Context, ToCompute};

pub fn fixpoint_only<
    Opt: Sync,
    const F: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    Context {
        fun,
        boundary,
        local_interaction,
        vs,
        k,
        integrated,
        r: (a, b),
        dt: dto,
        dx,
        maxdt,
        er,
        t,
        tend: _,
        opt,
    }: &mut Context<Opt, F, VX, VY, S>,
) -> f64 {
    *dto = maxdt.min(*dto);
    let [sizex, sizey] = *local_interaction;
    let mut err = 1.0;
    let mut cost = 0.0;
    let ko = *k;
    let mut fu = [[[[0.0f64; F]; VX]; VY]; S];
    let mut vdtk = [[[[0.0f64; F]; VX]; VY]; S];
    let mut errs = [[true; VX]; VY];
    let mut dt = *dto;
    let mut iter = 0;
    let maxiter = 10;
    let muldt: f64 = 1.01;
    let divdt: f64 = 0.5;
    let mut reset = false;
    while err > *er {
        iter += 1;
        if iter > maxiter {
            if reset {
                dt *= divdt;
            } else {
                dt /= muldt.powf(3.0);
            }
            reset = true;
            iter = 0;
            *k = ko;
            errs = [[true; VX]; VY];
        }
        for vy in 0..VY {
            for vx in 0..VX {
                for f in 0..F {
                    for s in 0..S {
                        if integrated[f] {
                            vdtk[s][vy][vx][f] = vs[vy][vx][f];
                            for s1 in 0..S {
                                vdtk[s][vy][vx][f] += dt * a[s][s1] * k[s1][vy][vx][f];
                            }
                        } else {
                            vdtk[s][vy][vx][f] = k[s][vy][vx][f];
                        }
                    }
                }
            }
        }

        err = 0.0;
        let ot = *t;
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let t = ot + cdt;
            fu[s]
                .par_iter_mut()
                .enumerate()
                .flat_map(|(vy, fsy)| {
                    fsy.par_iter_mut()
                        .enumerate()
                        .map(move |(vx, fsyx)| (vy, vx, fsyx))
                })
                .for_each(|(vy, vx, fu)| {
                    if errs[vy][vx] {
                        *fu = fun(
                            [&vs, &vdtk[s]],
                            boundary,
                            [vx as i32, vy as i32],
                            *dx,
                            *er,
                            [ot, t],
                            [dt, cdt],
                            opt,
                            ToCompute::All,
                        );
                    }
                });
        }
        for vy in 0..VY {
            for vx in 0..VX {
                if errs[vy][vx] {
                    cost += S as f64;
                    errs[vy][vx] = false;
                    for s in 0..S {
                        for f in 0..F {
                            let e = (fu[s][vy][vx][f] - k[s][vy][vx][f]).abs();
                            err = err.max(e);
                            errs[vy][vx] |= e > *er;
                        }
                    }
                }
            }
        }
        *k = fu;
        let mut tmperrs = [[false; VX]; VY];
        for vy in 0..VY {
            for vx in 0..VX {
                for dy in -sizey..=sizey {
                    for dx in -sizex..=sizex {
                        tmperrs[vy][vx] |=
                            errs[boundary[1](vy as i32 + dy, VY)][boundary[0](vx as i32 + dx, VX)];
                    }
                }
            }
        }
        errs = tmperrs;
    }
    let b = if let Some(b) = b { *b } else { a[S - 1] };
    for f in 0..F {
        for vy in 0..VY {
            for vx in 0..VX {
                if integrated[f] {
                    for s in 0..S {
                        vs[vy][vx][f] += dt * b[s] * k[s][vy][vx][f];
                    }
                } else {
                    vs[vy][vx][f] = k[S - 1][vy][vx][f];
                }
            }
        }
    }
    *t += dt;
    *dto = maxdt.min(dt * muldt);
    cost / (VX * VY) as f64
}