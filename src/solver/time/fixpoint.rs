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
        er,
        t,
        ot,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
        freezeout_energy: _,
    }: &mut Context<Opt, F, C, VX, VY, S>,
    err_ref_p: Option<(usize, f64)>,
) -> Option<(f64, [[usize; VX]; VY], usize)> {
    let [sizex, sizey] = *local_interaction;
    *dto = maxdt.min(*dto);
    let mut iserr = true;
    let ko = *k;
    let mut fu = *k;
    let mut vdtk = [[[0.0f64; F]; VX]; VY];
    let mut trdtk = [[[0.0f64; C]; VX]; VY];
    let mut errs = [[true; VX]; VY];
    let mut nbiter = [[0usize; VX]; VY];
    let mut dt = *dto;
    let mut iter = 0;
    let mut failed = 1;
    let maxiter = 10;
    let maxfailed = 30;
    let divdt: f64 = 0.5;
    while iserr {
        iter += 1;
        if iter > failed * maxiter {
            eprintln!("fail {}", failed);
            failed += 1;
            if failed >= maxfailed {
                return None;
            }
            dt *= divdt;
            iter = 0;
            *k = ko;
            errs = [[true; VX]; VY];
        }
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let ct = *t + cdt;
            pfor2d2(&mut vdtk, &mut trdtk, &|(Coord { x, y }, vdtk, trdtk)| {
                for f in 0..F {
                    vdtk[f] = vs[y][x][f];
                    for s1 in 0..S {
                        vdtk[f] += dt * a[s][s1] * fu[s1][y][x][f];
                    }
                    (*vdtk, *trdtk) = constraints(ct, *vdtk, *trdtk);
                }
            });
            pfor2d(&mut fu[s], &|(Coord { x, y }, fu)| {
                if errs[y][x] {
                    *fu = fun(
                        [&vs, &vdtk],
                        [&trs, &trdtk],
                        constraints,
                        boundary,
                        [x as i32, y as i32],
                        *dx,
                        [*t, ct],
                        [dt, cdt],
                        opt,
                    );
                }
            });
        }

        let err_ref: Box<dyn Fn(usize, usize) -> f64 + Sync> = if let Some((_f, mi)) = err_ref_p {
            // let trs = &trs;
            // let m = trs.iter().fold(0.0, |acc: f64, v| {
            //     acc.max(v.iter().fold(0.0, |acc, v| acc.max(v[f])))
            // });
            // Box::new(move |x, y| (trs[y][x][f] / m).powi(2)) // due to '/m', the reference is in <=1, WARNING taking the square might be problematic
            // Box::new(move |x, y| (trs[y][x][f] / m).powi(2).max(mi))
            Box::new(move |_, _| mi)
        } else {
            Box::new(|_, _| 1.0)
        };

        pfor2d2(&mut errs, &mut nbiter, &|(Coord { x, y }, errs, nbiter)| {
            if *errs {
                *errs = false;
                *nbiter += 1;
                for s in 0..S {
                    for f in 0..F {
                        let e = (fu[s][y][x][f] - k[s][y][x][f]).abs();
                        if e.is_nan() {
                            panic!("NaN");
                        }
                        *errs |= e * err_ref(x, y) > *er || e.is_nan(); // |f(k)-k| * cutoff > dt^p
                                                                        // *errs |= e * err_ref(x, y)
                                                                        //     > (wb::T(trs[y][x][f]).max(20.0 / 200.0).powi(2) * *er)
                                                                        //     || e.is_nan();
                                                                        // |f(k)-k| * cutoff > dt^p
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
    *ovs = *vs;
    *otrs = *trs;
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
    *ot = *t;
    *t += dt;
    *dto = maxdt.min(dt * 1.1);
    for vy in 0..VY {
        for vx in 0..VX {
            let tmp = vs[vy][vx];
            (vs[vy][vx], trs[vy][vx]) = constraints(*t, vs[vy][vx], trs[vy][vx]);
            for f in 0..F {
                total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
            }
        }
    }
    if let Some(post) = post_constraints {
        for vy in 0..VY {
            for vx in 0..VX {
                let tmp = vs[vy][vx];
                (vs[vy][vx], trs[vy][vx]) = post(*t, vs[vy][vx], trs[vy][vx]);
                for f in 0..F {
                    total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
                }
            }
        }
    }
    Some((cost / (VX * VY) as f64, nbiter, failed - 1))
}
