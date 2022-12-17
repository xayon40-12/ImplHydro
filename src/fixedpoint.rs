use rayon::prelude::*;

use crate::newton::Context;

pub fn fixedpoint<Opt: Sync, const F: usize, const VX: usize, const VY: usize, const S: usize>(
    Context {
        fun,
        boundary,
        local_interaction,
        vs,
        k,
        integrated,
        r,
        dt: dto,
        dx,
        maxdt,
        er,
        t,
        tend: _,
        opt,
    }: &mut Context<Opt, F, VX, VY, S>,
) -> usize {
    *dto = maxdt.min(*dto);
    let [_sizex, _sizey] = *local_interaction;
    let mut err = 1.0;
    let mut iterations = 0;
    let ko = *k;
    let mut fu = [[[[0.0f64; F]; VX]; VY]; S];
    let mut vdtk = [[[[0.0f64; F]; VX]; VY]; S];
    let mut errs = [[[true; VY]; VY]; S];
    let mut dt = *dto;
    let mut iter = 0;
    let maxiter = 10;
    while err > *er {
        iterations += 1;
        iter += 1;
        if iter > maxiter {
            dt *= 0.5;
            iter = 0;
            *k = ko;
            errs = [[[true; VY]; VY]; S];
        }
        for f in 0..F {
            for vy in 0..VY {
                for vx in 0..VX {
                    for s in 0..S {
                        if integrated[f] {
                            vdtk[s][vy][vx][f] = 0.0;
                            for s1 in 0..S {
                                vdtk[s][vy][vx][f] +=
                                    vs[vy][vx][f] + dt * r[s][s1] * k[s1][vy][vx][f];
                            }
                        } else {
                            vdtk[s][vy][vx][f] = k[s][vy][vx][f];
                        }
                    }
                }
            }
        }

        err = 0.0;
        for s in 0..S {
            let ot = *t;
            let c = r[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let t = ot + cdt;
            for vy in 0..VY {
                fu[s][vy].par_iter_mut().enumerate().for_each(|(vx, fu)| {
                    *fu = fun(
                        [&vs, &vdtk[s]],
                        boundary,
                        [vx as i32, vy as i32],
                        *dx,
                        [ot, t],
                        [dt, cdt],
                        opt,
                    );
                });
                for vx in 0..VX {
                    for f in 0..F {
                        let e = (fu[s][vy][vx][f] - k[s][vy][vx][f]).abs();
                        err = err.max(e);
                        errs[s][vy][vx] = e <= *er;
                    }
                    k[s][vy][vx] = fu[s][vy][vx];
                }
            }
        }
        // eprintln!("i: {}, err: {}", iterations, err);
    }
    for f in 0..F {
        for vy in 0..VY {
            for vx in 0..VX {
                if integrated[f] {
                    for s in 0..S {
                        vs[vy][vx][f] += dt * r[S - 1][s] * k[s][vy][vx][f];
                    }
                } else {
                    vs[vy][vx][f] = k[S - 1][vy][vx][f];
                }
            }
        }
    }
    *t += dt;
    *dto = maxdt.min(dt * 1.1);
    iterations
}
