use rayon::prelude::*;

use crate::solver::{
    context::Context,
    utils::{pfor2d, pfor2d2, Coord},
};

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
        post_constraints,
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
    let mut iserr = true;
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
    while iserr {
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
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let t = *ot + cdt;
            pfor2d(&mut vdtk, &|(Coord { x, y }, vdtk)| {
                for f in 0..F {
                    vdtk[f] = vs[y][x][f];
                    for s1 in 0..S {
                        vdtk[f] += dt * a[s][s1] * fu[s1][y][x][f];
                    }
                    *vdtk = constraints(t, *vdtk);
                }
            });
            pfor2d(&mut fu[s], &|(Coord { x, y }, fu)| {
                if errs[y][x] {
                    *fu = fun(
                        [&vs, &vdtk],
                        constraints,
                        transform,
                        boundary,
                        [x as i32, y as i32],
                        *dx,
                        [*ot, t],
                        [dt, cdt],
                        opt,
                    );
                }
            });
        }

        pfor2d2(&mut errs, &mut nbiter, &|(Coord { x, y }, errs, nbiter)| {
            if *errs {
                *errs = false;
                *nbiter += 1;
                for s in 0..S {
                    for f in 0..F {
                        let e = (fu[s][y][x][f] - k[s][y][x][f]).abs();
                        *errs |= e > *er;
                    }
                }
            }
        });

        iserr = errs
            .par_iter()
            .flat_map(|e| e.par_iter())
            .map(|v| *v)
            .reduce(|| false, |acc, a| acc || a);
        *k = fu;
        let mut tmperrs = [[false; VX]; VY];
        pfor2d(&mut tmperrs, &|(Coord { x, y }, tmperrs)| {
            for dy in -sizey..=sizey {
                for dx in -sizex..=sizex {
                    *tmperrs |=
                        errs[boundary[1](y as i32 + dy, VY)][boundary[0](x as i32 + dx, VX)];
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
    let cost = (S * nbiter
        .par_iter()
        .flat_map(|e| e.par_iter())
        .map(|v| *v)
        .reduce(|| 0, |acc, a| acc + a)) as f64;
    *ot += dt;
    *dto = maxdt.min(dt * 1.1);
    if let Some(post) = post_constraints {
        for vy in 0..VY {
            for vx in 0..VX {
                vs[vy][vx] = post(*ot, vs[vy][vx]);
            }
        }
    }
    Some((cost / (VX * VY) as f64, nbiter))
}
