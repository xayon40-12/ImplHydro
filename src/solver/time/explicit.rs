use crate::solver::{context::Context, pfor2d};

use super::schemes::Scheme;

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
        transform,
        boundary,
        local_interaction: _,
        vs,
        k, // k[S-1] is used as old vs
        r: Scheme { aij: a, .. },
        dt,
        dx,
        maxdt,
        er: _,
        t: ot,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
    }: &mut Context<Opt, F, C, VX, VY, S>,
) -> Option<(f64, [[usize; VX]; VY])> {
    *dt = maxdt.min(*dt);
    let cost = S as f64;
    let nbiter = [[1usize; VX]; VY];
    let mut fu = [[[0.0f64; F]; VX]; VY];
    let mut vdtk = *vs;

    let mut c;
    let mut cdt = *dt;
    let mut t = *ot;
    for s in 0..S {
        for vy in 0..VY {
            for vx in 0..VX {
                vdtk[vy][vx] = constraints(t, vdtk[vy][vx]);
            }
        }
        pfor2d(&mut fu, &|(vy, vx, fu)| {
            *fu = fun(
                [&k[S - 1], &vdtk],
                constraints,
                transform,
                boundary,
                [vx as i32, vy as i32],
                *dx,
                [*ot, t],
                [*dt, cdt],
                opt,
            );
        });
        k[S - 1] = vdtk;
        k[s] = fu;
        pfor2d(&mut vdtk, &|(vy, vx, vdtk)| {
            for f in 0..F {
                vdtk[f] = vs[vy][vx][f];
                for s1 in 0..S {
                    vdtk[f] += *dt * a[s][s1] * k[s1][vy][vx][f];
                }
            }
        });
        c = a[s].iter().fold(0.0, |acc, r| acc + r);
        cdt = t;
        t = *ot + c * *dt;
        cdt = t - cdt; // difference between current time and previous time
    }
    k[S - 1] = *vs; // store the old vs in k[S-1] for next time step
    *vs = vdtk;

    *ot += *dt;
    Some((cost, nbiter))
}
