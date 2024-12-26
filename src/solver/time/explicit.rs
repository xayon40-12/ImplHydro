use crate::{
    solver::{
        context::Context,
        utils::{cfor3d, gen_coords, Coord},
    },
    FLOAT,
};
use boxarray::boxarray;

use super::schemes::Scheme;

pub fn explicit<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
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
        dxs,
        maxdt,
        ot,
        t,
        t0: _,
        tend: _,
        opt,
        p: _,
        dpde: _,
        freezeout_energy: _,
    }: &mut Context<Opt, F, C, VX, VY, VZ, S>,
) -> Option<(FLOAT, Box<[[[usize; VX]; VY]; VZ]>)> {
    let coords = gen_coords::<VX, VY, VZ>();
    *dt = maxdt.min(*dt);
    let cost = S as FLOAT;
    let nbiter: Box<[[[usize; VX]; VY]; VZ]> = boxarray(1);
    let mut fu: Box<[[[[FLOAT; F]; VX]; VY]; VZ]> = boxarray(0.0);
    let mut vdtk = vs.clone();
    let mut trdtk = trs.clone();
    let mut ovdtk = ovs.clone();
    let mut otrdtk = otrs.clone();
    let mut ct = *t;

    let mut c;
    let mut cdt = *t - *ot;
    for s in 0..S {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    (vdtk[vz][vy][vx], trdtk[vz][vy][vx]) = constraints(ct, vdtk[vz][vy][vx]);
                }
            }
        }
        let es = if s == 0 { S - 1 } else { s - 1 }; // index of last computed k
        cfor3d(&coords, &mut fu, |&Coord { x, y, z, .. }, fu| {
            *fu = fun(
                &k[es],
                [&ovdtk, &vdtk],
                [&otrdtk, &trdtk],
                constraints,
                boundary,
                [x as i32, y as i32, z as i32],
                *dxs,
                [*ot, ct],
                [*dt, cdt],
                opt,
            );
        });
        otrdtk = trdtk.clone();
        ovdtk = vdtk.clone();
        k[s] = *fu;
        cfor3d(&coords, &mut vdtk, |&Coord { x, y, z, .. }, vdtk| {
            for f in 0..F {
                vdtk[f] = vs[z][y][x][f];
                for s1 in 0..S {
                    vdtk[f] += *dt * a[s][s1] * k[s1][z][y][x][f];
                }
            }
        });
        c = a[s].iter().fold(0.0, |acc, r| acc + r);
        cdt = ct;
        *ot = ct;
        ct = *t + c * *dt;
        cdt = ct - cdt; // difference between current time and previous time
    }
    *ovs = vs.clone();
    *otrs = trs.clone();
    *vs = vdtk;

    *ot = *t;
    *t += *dt;
    for vz in 0..VZ {
        for vy in 0..VY {
            for vx in 0..VX {
                let tmp = vs[vz][vy][vx];
                (vs[vz][vy][vx], trs[vz][vy][vx]) = constraints(*t, vs[vz][vy][vx]);
                for f in 0..F {
                    total_diff_vs[vz][vy][vx][f] += (tmp[f] - vs[vz][vy][vx][f]).abs();
                }
            }
        }
    }
    if let Some(post) = post_constraints {
        for vz in 0..VZ {
            for vy in 0..VY {
                for vx in 0..VX {
                    let tmp = vs[vz][vy][vx];
                    (vs[vz][vy][vx], trs[vz][vy][vx]) = post(*t, vs[vz][vy][vx]);
                    for f in 0..F {
                        total_diff_vs[vz][vy][vx][f] += (tmp[f] - vs[vz][vy][vx][f]).abs();
                    }
                }
            }
        }
    }
    Some((cost, nbiter))
}
