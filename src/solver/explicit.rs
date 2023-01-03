use rayon::prelude::*;

use crate::solver::context::{Context, ToCompute};

pub fn explicit<
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
        boundary,
        local_interaction: _,
        vs,
        k,
        integrated,
        r: (a, _),
        dt,
        dx,
        maxdt,
        er,
        t: ot,
        t0: _,
        tend: _,
        opt,
    }: &mut Context<Opt, F, C, VX, VY, S>,
) -> f64 {
    *dt = maxdt.min(*dt);
    let cost = S as f64;
    let mut fu = [[[0.0f64; F]; VX]; VY];
    let mut vdtk = *vs;

    let mut c;
    let mut cdt = 0.0;
    let mut t = *ot;
    for s in 0..S {
        fu.par_iter_mut()
            .enumerate()
            .flat_map(|(vy, fsy)| {
                fsy.par_iter_mut()
                    .enumerate()
                    .map(move |(vx, fsyx)| (vy, vx, fsyx))
            })
            .for_each(|(vy, vx, fu)| {
                *fu = fun(
                    [&vs, &vdtk],
                    constraints,
                    boundary,
                    [vx as i32, vy as i32],
                    *dx,
                    *er,
                    [*ot, t],
                    [*dt, cdt],
                    opt,
                    ToCompute::Integrated,
                );
            });
        for vy in 0..VY {
            for vx in 0..VX {
                for f in 0..F {
                    if integrated[f] {
                        vdtk[vy][vx][f] = vs[vy][vx][f];
                        k[s][vy][vx][f] = fu[vy][vx][f];
                        for s1 in 0..S {
                            vdtk[vy][vx][f] += *dt * a[s][s1] * k[s1][vy][vx][f];
                        }
                    }
                }
            }
        }
        c = a[s].iter().fold(0.0, |acc, r| acc + r);
        cdt = c * *dt;
        t = *ot + cdt;
        fu.par_iter_mut()
            .enumerate()
            .flat_map(|(vy, fsy)| {
                fsy.par_iter_mut()
                    .enumerate()
                    .map(move |(vx, fsyx)| (vy, vx, fsyx))
            })
            .for_each(|(vy, vx, fu)| {
                *fu = fun(
                    [&vs, &vdtk],
                    constraints,
                    boundary,
                    [vx as i32, vy as i32],
                    *dx,
                    *er,
                    [*ot, t],
                    [*dt, cdt],
                    opt,
                    ToCompute::NonIntegrated,
                );
            });
        for vy in 0..VY {
            for vx in 0..VX {
                for f in 0..F {
                    if !integrated[f] {
                        k[s][vy][vx][f] = fu[vy][vx][f];
                        vdtk[vy][vx][f] = k[s][vy][vx][f];
                    }
                }
            }
        }
    }

    let b = a[S - 1];
    for f in 0..F {
        for vy in 0..VY {
            for vx in 0..VX {
                if integrated[f] {
                    for s in 0..S {
                        vs[vy][vx][f] += *dt * b[s] * k[s][vy][vx][f];
                    }
                } else {
                    vs[vy][vx][f] = k[S - 1][vy][vx][f];
                }
            }
        }
    }
    *ot += *dt;
    cost
}
