use rayon::prelude::*;

use crate::solver::{
    context::Context,
    utils::{pfor2d, pfor2d2, Coord},
};

use super::schemes::Scheme;

#[derive(Clone)]
pub struct FixCTX {
    pub fka: f64,
}

pub fn fixpoint<
    Opt: Clone + Sync,
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
    err_c: Option<f64>,
    ctx: &mut FixCTX,
) -> Option<(f64, [[usize; VX]; VY], usize)> {
    let err_c = err_c.unwrap_or(1.0);
    let [sizex, sizey] = *local_interaction;
    *dto = maxdt.min(*dto);
    let mut iserr = true;
    let ko = *k;
    let mut tmpk = *k;
    let mut fu = *k;
    let mut tmpfu = *k;
    let mut vdtk = [[[0.0f64; F]; VX]; VY];
    let mut trdtk = [[[0.0f64; C]; VX]; VY];
    let mut errs = [[true; VX]; VY];
    let mut nbiter = [[0usize; VX]; VY];
    let dt = *dto;
    let mut iter = 0;
    let mut failed = 1usize;
    let maxiter = 10;
    let maxfailed = 30;
    let mut maxe = 1e50f64;
    let mut omaxe = maxe;
    // let mut reset_ffka = ctx.fka;
    let mut reset_ffka = 1.0;
    // let fkafact = 0.99f64;
    while iserr {
        iter += 1;
        // let mul = 1 + (reset_ffka.ln() / fkafact.ln()) as usize;
        // println!("mul: {}, fka: {:.3e}", mul, reset_ffka);
        let mul = failed;
        if iter > mul * maxiter || maxe > 1e50 {
            eprintln!("fail {}", failed);
            failed += 1;
            // if reset_ffka > fkafact {
            //     reset_ffka = fkafact;
            // }
            // reset_ffka *= reset_ffka;
            reset_ffka *= 0.9;
            maxe = 1e50;
            omaxe = maxe;
            iter = 0;
            *k = ko;
            fu = *k;
            tmpfu = fu;
            vdtk = [[[0.0f64; F]; VX]; VY];
            trdtk = [[[0.0f64; C]; VX]; VY];
            errs = [[true; VX]; VY];
            if failed >= maxfailed {
                return None;
            }
        }
        for s in 0..S {
            let c = a[s].iter().fold(0.0, |acc, r| acc + r);
            let cdt = c * dt;
            let ct = *t + cdt;
            pfor2d2(&mut vdtk, &mut trdtk, &|(Coord { x, y }, vdtk, trdtk)| {
                for f in 0..F {
                    vdtk[f] = vs[y][x][f];
                    for s1 in 0..S {
                        vdtk[f] += dt * a[s][s1] * k[s1][y][x][f];
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

        pfor2d2(&mut errs, &mut nbiter, &|(Coord { x, y }, errs, nbiter)| {
            if *errs {
                *errs = false;
                *nbiter += 1;
                for s in 0..S {
                    for f in 0..F {
                        let fuf = fu[s][y][x][f];
                        let kf = k[s][y][x][f];
                        let e = (fuf - kf).abs();
                        if e.is_nan() {
                            panic!(
                                "\n\nNaN encountered in fixpoint iteration.\nfields: {:?}\n\n",
                                fu[s][y][x]
                            );
                        }
                        *errs |= e > err_c * *er || e.is_nan(); // |f(k)-k| * cutoff > dt^p
                                                                // *errs |= e > 1e-2 || e.is_nan();
                    }
                }
            }
        });

        iserr = errs
            .par_iter()
            .flat_map(|e| e.par_iter())
            .map(|v| *v)
            .reduce(|| false, |acc, a| acc || a);

        maxe = 0.0;
        for s in 0..S {
            for j in 0..VY {
                for i in 0..VX {
                    for f in 0..F {
                        let kk = k[s][j][i][f];
                        let ff = fu[s][j][i][f];
                        let m = (kk - ff).abs();
                        maxe = maxe.max(m);
                    }
                }
            }
        }

        let enable_fka = true;
        let mut ffka = reset_ffka;
        if maxe >= omaxe && enable_fka {
            *k = tmpk;
            fu = tmpfu;
            ffka = 0.9 * omaxe / maxe;
        } else {
            tmpfu = fu;
            tmpk = *k;
            omaxe = maxe;
        }
        for s in 0..S {
            for j in 0..VY {
                for i in 0..VX {
                    for f in 0..F {
                        let kk = k[s][j][i][f];
                        let ff = fu[s][j][i][f];
                        k[s][j][i][f] = kk + ffka * (ff - kk);
                    }
                }
            }
        }
        let debug = false;
        // let debug = true;
        if debug {
            println!(
                "t: {:.3e}, ffka: {:.3e}, max: {:.3e}, omax: {:.3e}, lim: {:.3e}",
                t,
                ffka,
                maxe,
                omaxe,
                err_c * *er
            );
        }

        let mut tmperrs = [[false; VX]; VY];
        pfor2d(&mut tmperrs, &|(Coord { x, y }, tmperrs)| {
            let mut update = |dx, dy| {
                *tmperrs |= errs[boundary[1](y as i32 + dy, VY)][boundary[0](x as i32 + dx, VX)]
            };
            for dy in -sizey..=sizey {
                if dy == 0 {
                    for dx in -sizex..=sizex {
                        update(dx, dy);
                    }
                } else {
                    update(0, dy);
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
    ctx.fka = if failed == 1 {
        reset_ffka.sqrt()
    } else {
        reset_ffka
    };
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
