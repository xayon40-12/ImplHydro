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

pub type ErrThr<'a, const F: usize, const C: usize, const VX: usize, const VY: usize> =
    &'a dyn Fn(f64, &[[[f64; F]; VX]; VY], &[[[f64; C]; VX]; VY]) -> f64;

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
        t,
        ot,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
        freezeout_energy: _,
    }: &mut Context<Opt, F, C, VX, VY, S>,
    err_thr: ErrThr<F, C, VX, VY>,
    _ctx: &mut FixCTX,
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
    let dt = *dto;
    let mut failed = 1usize;
    let maxfailed = 30;
    let mut maxe = 1e50f64;
    let mut omaxe = maxe;
    // let mut ffka = ctx.fka;
    let mut ffka = 1.0;
    let mut maxe_count = 0;
    while iserr {
        if maxe > 1e50 || maxe_count > 4 {
            eprintln!("fail {}", failed);
            failed += 1;
            // ffka *= 0.9;
            ffka *= (*dx / dt).min(0.5);
            maxe = 1e50;
            omaxe = maxe;
            *k = ko;
            fu = *k;
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
                    (*vdtk, *trdtk) = constraints(ct, *vdtk);
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

        let er = err_thr(*t, &vs, &trs);
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
                        *errs |= e > er || e.is_nan(); // |f(k)-k| * cutoff > dt^p
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

        // omaxe = maxe;
        maxe = 0.0;
        let mut errr = [[0.0f64; VX]; VY];
        for s in 0..S {
            for j in 0..VY {
                for i in 0..VX {
                    for f in 0..F {
                        let kk = k[s][j][i][f];
                        let ff = fu[s][j][i][f];
                        let d = ff - kk;
                        errr[j][i] = errr[j][i].max(d.abs());
                        maxe = maxe.max(d.abs());
                        k[s][j][i][f] = kk + ffka * d;
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
            let miter = nbiter.iter().flat_map(|v| v.iter()).max().unwrap();
            println!(
                "t: {:.3e}, ffka: {:.3e}, maxe: {:.3e}, omaxe: {:.3e}, lim: {:.3e}, miter: {}",
                t, ffka, maxe, omaxe, er, miter
            );

            let col: Vec<char> = " .:-=+*%#@".chars().collect();
            let mc = col.len() - 1;
            let iter = nbiter
                .iter()
                .map(|v| {
                    v.iter()
                        .map(|v| format!("{s}{s}", s = col[v * mc / miter].to_string()))
                        .collect::<Vec<_>>()
                        .join("")
                })
                .collect::<Vec<_>>();
            let err = errr
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

        let mut tmperrs = [[false; VX]; VY];
        pfor2d(&mut tmperrs, &|(Coord { x, y }, tmperrs)| {
            let mut update = |dx, dy| {
                let i = x as i32 + dx;
                let j = y as i32 + dy;
                if i >= 0 && i < VX as i32 && j >= 0 && j < VY as i32 {
                    *tmperrs |= errs[j as usize][i as usize];
                }
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
    // println!("\nnext");
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
    // ctx.fka = if failed == 1 { ffka.sqrt() } else { ffka };
    for vy in 0..VY {
        for vx in 0..VX {
            let tmp = vs[vy][vx];
            (vs[vy][vx], trs[vy][vx]) = constraints(*t, vs[vy][vx]);
            for f in 0..F {
                total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
            }
        }
    }
    if let Some(post) = post_constraints {
        for vy in 0..VY {
            for vx in 0..VX {
                let tmp = vs[vy][vx];
                (vs[vy][vx], trs[vy][vx]) = post(*t, vs[vy][vx]);
                for f in 0..F {
                    total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
                }
            }
        }
    }
    Some((cost / (VX * VY) as f64, nbiter, failed - 1))
}
