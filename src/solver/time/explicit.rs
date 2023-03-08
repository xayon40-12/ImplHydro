use crate::solver::{
    context::Context,
    utils::{pfor2d, Coord},
};

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
        boundary,
        post_constraints,
        local_interaction: _,
        total_diff_vs,
        vstrs: (vs, trs),
        ovstrs: (ovs, otrs),
        k,
        r: Scheme { aij: a, .. },
        dt,
        dx,
        maxdt,
        er: _,
        ot,
        t,
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
    let mut trdtk = *trs;
    let mut ovdtk = *ovs;
    let mut otrdtk = *otrs;
    let mut ct = *t;

    let mut c;
    let mut cdt = *t - *ot;
    for s in 0..S {
        for vy in 0..VY {
            for vx in 0..VX {
                (vdtk[vy][vx], trdtk[vy][vx]) = constraints(ct, vdtk[vy][vx], trdtk[vy][vx]);
            }
        }
        pfor2d(&mut fu, &|(Coord { x, y }, fu)| {
            *fu = fun(
                [&ovdtk, &vdtk],
                [&otrdtk, &trdtk],
                constraints,
                boundary,
                [x as i32, y as i32],
                *dx,
                [*ot, ct],
                [*dt, cdt],
                opt,
            );
        });
        otrdtk = trdtk;
        ovdtk = vdtk;
        k[s] = fu;
        pfor2d(&mut vdtk, &|(Coord { x, y }, vdtk)| {
            for f in 0..F {
                vdtk[f] = vs[y][x][f];
                for s1 in 0..S {
                    vdtk[f] += *dt * a[s][s1] * k[s1][y][x][f];
                }
            }
        });
        c = a[s].iter().fold(0.0, |acc, r| acc + r);
        cdt = ct;
        *ot = ct;
        ct = *t + c * *dt;
        cdt = ct - cdt; // difference between current time and previous time
    }
    *ovs = *vs; // store the old vs in k[S-1] for next time step
    *otrs = *trs;
    *vs = vdtk;

    *ot = *t;
    *t += *dt;
    for vy in 0..VY {
        for vx in 0..VX {
            let tmp = vs[vy][vx];
            (vs[vy][vx], trs[vy][vx]) = constraints(*ot, vs[vy][vx], trs[vy][vx]);
            for f in 0..F {
                total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
            }
        }
    }
    if let Some(post) = post_constraints {
        for vy in 0..VY {
            for vx in 0..VX {
                let tmp = vs[vy][vx];
                (vs[vy][vx], trs[vy][vx]) = post(*ot, vs[vy][vx], trs[vy][vx]);
                for f in 0..F {
                    total_diff_vs[vy][vx][f] += (tmp[f] - vs[vy][vx][f]).abs();
                }
            }
        }
    }
    Some((cost, nbiter))
}
