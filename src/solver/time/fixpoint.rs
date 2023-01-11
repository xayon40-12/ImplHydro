use crate::solver::{context::Context, pfor2d};

use super::schemes::Scheme;

pub fn fixpoint<
    Opt: Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    Context {
        fun,
        constraints,
        transform,
        boundary,
        local_interaction,
        vs,
        k,
        r: Scheme { aij: a, bj: b, .. },
        dt: dto,
        dx,
        maxdt,
        er,
        t: ot,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
    }: &mut Context<Opt, F, C, VX, VY, S>,
) -> Option<(f64, [[usize; VX]; VY])> {
    let [sizex, sizey] = *local_interaction;
    *dto = maxdt.min(*dto);
    let mut err = 1.0;
    let mut cost = 0.0;
    let ko = *k;
    let mut fu = *k;
    let mut vdtk = [[[0.0f64; F]; VX]; VY];
    let mut errs = [[true; VX]; VY];
    let mut nbiter = [[0usize; VX]; VY];
    let mut dt = *dto;
    let mut iter = 0;
    let mut failed = 0;
    let maxiter = 50;
    let maxfailed = 10;
    let muldt: f64 = 1.01;
    let divdt: f64 = 0.5;
    let mut reset = false;
    while err > *er {
        iter += 1;
        if iter > maxiter {
            failed += 1;
            if failed >= maxfailed {
                return None;
            }
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
        err = 0.0;
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let t = *ot + cdt;
            for vy in 0..VY {
                for vx in 0..VX {
                    for f in 0..F {
                        vdtk[vy][vx][f] = vs[vy][vx][f];
                        for s1 in 0..S {
                            vdtk[vy][vx][f] += dt * a[s][s1] * fu[s1][vy][vx][f];
                        }
                        vdtk[vy][vx] = constraints(t, vdtk[vy][vx]);
                    }
                }
            }
            pfor2d(&mut fu[s], &|(vy, vx, fu)| {
                if errs[vy][vx] {
                    *fu = fun(
                        [&vs, &vdtk],
                        constraints,
                        transform,
                        boundary,
                        [vx as i32, vy as i32],
                        *dx,
                        [*ot, t],
                        [dt, cdt],
                        opt,
                    );
                }
            });
        }

        for vy in 0..VY {
            for vx in 0..VX {
                if errs[vy][vx] {
                    cost += S as f64;
                    errs[vy][vx] = false;
                    nbiter[vy][vx] += 1;
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
        pfor2d(&mut tmperrs, &|(vy, vx, tmperrs)| {
            for dy in -sizey..=sizey {
                for dx in -sizex..=sizex {
                    *tmperrs |=
                        errs[boundary[1](vy as i32 + dy, VY)][boundary[0](vx as i32 + dx, VX)];
                }
            }
        });
        errs = tmperrs;
    }
    let b = if let Some(b) = b { *b } else { a[S - 1] };
    for f in 0..F {
        for vy in 0..VY {
            for vx in 0..VX {
                for s in 0..S {
                    vs[vy][vx][f] += dt * b[s] * k[s][vy][vx][f];
                }
            }
        }
    }
    *ot += dt;
    *dto = maxdt.min(dt * 1.1);
    Some((cost / (VX * VY) as f64, nbiter))
}
