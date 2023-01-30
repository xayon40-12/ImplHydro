use crate::solver::{context::Boundary, Transform};

use super::{Dir, Eigenvalues, Flux};

pub fn kl<const F: usize, const VX: usize, const VY: usize, const C: usize>(
    vars: &[[[f64; F]; VX]; VY],
    [boundx, boundy]: &[Boundary; 2],
    [x, y]: [i32; 2],
    dir: Dir,
    t: f64,
    flux: [Flux<C>; F],
    constraints: Transform<F, F>,
    transform: Transform<F, C>,
    eigenvalues: Eigenvalues<C>,
    pre_flux_limiter: Transform<F, F>,
    post_flux_limiter: Transform<F, F>,
    dx: f64,
    p: f64, // in kl 'p' is the equivalent of 'theta' in the minmod flux limiter of kt
) -> [f64; F] {
    let mut vals: [[f64; F]; 5] = [[0.0; F]; 5];
    let mut tvals: [[f64; C]; 5] = [[0.0; C]; 5];
    for l in 0..5 {
        let (lx, ly) = match dir {
            Dir::X => (l as i32 - 2, 0),
            Dir::Y => (0, l as i32 - 2),
        };
        vals[l] = constraints(t, vars[boundy(y + ly, VY)][boundx(x + lx, VX)]);
        tvals[l] = transform(t, vals[l]);
        vals[l] = pre_flux_limiter(t, vals[l]);
    }
    let mut valso: [[[f64; F]; 3]; 3] = [[[0.0; F]; 3]; 3];
    for l in 0..3 {
        let (lx, ly) = match dir {
            Dir::X => (l as i32 - 1, 0),
            Dir::Y => (0, l as i32 - 1),
        };
        for o in 0..3 {
            let (ox, oy) = match dir {
                Dir::X => (0, l as i32 - 1),
                Dir::Y => (l as i32 - 1, 0),
            };
            valso[l][o] = pre_flux_limiter(
                t,
                constraints(t, vars[boundy(y + ly + oy, VY)][boundx(x + lx + ox, VX)]),
            );
        }
    }
    let eps: f64 = 1e-14;
    let mut is: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut diff: [[[f64; F]; 3]; 3] = [[[0.0; F]; 3]; 3];
    let mut diff2: [[f64; F]; 3] = [[0.0; F]; 3];
    let mut diff2o: [[f64; F]; 3] = [[0.0; F]; 3];
    let mut fj: [[f64; F]; 3] = [[0.0; F]; 3];
    for l in 0..3 {
        for f in 0..F {
            diff[l][0][f] = vals[l + 1][f] - vals[l][f];
            diff[l][1][f] = vals[l + 2][f] - vals[l][f];
            diff[l][2][f] = vals[l + 2][f] - vals[l + 1][f];

            diff2[l][f] = diff[l][2][f] - diff[l][0][f];
            diff2o[l][f] = valso[l][2][f] - 2.0 * valso[l][1][f] + valso[l][0][f]; // this is automatically 0 in 1D

            is[l][0] += diff[l][0][f].powi(2);
            is[l][1] += 13.0 / 3.0 * diff2[l][f].powi(2) + 1.0 / 4.0 * diff[l][1][f].powi(2);
            is[l][2] += diff[l][2][f].powi(2);

            fj[l][f] = flux[f](t, tvals[l + 1]);
        }
    }
    let c: [f64; 3] = [1.0 / 4.0, 1.0 / 2.0, 1.0 / 4.0];
    let mut alpha: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for l in 0..3 {
        for m in 0..3 {
            alpha[l][m] = c[m] / (eps + is[l][m]).powf(p);
        }
    }
    let mut sa: [f64; 3] = [0.0; 3];
    for l in 0..3 {
        for m in 0..3 {
            sa[l] += alpha[l][m];
        }
    }
    let mut w: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for l in 0..3 {
        for m in 0..3 {
            w[l][m] = alpha[l][m] / sa[l];
        }
    }

    let mut aa: [[f64; F]; 3] = [[0.0; F]; 3];
    let mut bb: [[f64; F]; 3] = [[0.0; F]; 3];
    let mut cc: [[f64; F]; 3] = [[0.0; F]; 3];
    for l in 0..3 {
        for f in 0..F {
            aa[l][f] = vals[l + 1][f] - w[l][1] / 12.0 * (diff2[l][f] + diff2o[l][f]);
            bb[l][f] =
                w[l][2] * diff[l][2][f] + w[l][1] * diff[l][1][f] / 2.0 + w[l][0] * diff[l][0][f]; // remove '/dx' as it is canceled in upm
            cc[l][f] = w[l][1] * diff2[l][f]; // remove '*2.0/dx2' as it is canceled in upm
        }
    }

    let mut upm: [[[f64; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
    let mut tupm: [[[f64; C]; 2]; 2] = [[[0.0; C]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            let s = pm as f64 - 0.5;
            let j = jpm + pm;
            for f in 0..F {
                upm[jpm][pm][f] = aa[j][f] - s * bb[j][f] + 0.25 * cc[j][f];
                // remove '*dx' for bb and '*dx2/2.0' for cc
            }
            upm[jpm][pm] = constraints(t, post_flux_limiter(t, upm[jpm][pm]));
            tupm[jpm][pm] = transform(t, upm[jpm][pm]);
        }
    }
    let mut fpm: [[[f64; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            for n in 0..F {
                fpm[jpm][pm][n] = flux[n](t, tupm[jpm][pm]);
            }
        }
    }
    let mut a: [f64; 2] = [0.0; 2];
    match eigenvalues {
        Eigenvalues::Analytical(eig_analytical) => {
            for jpm in 0..2 {
                a[jpm] = eig_analytical(t, tupm[jpm][0]).max(eig_analytical(t, tupm[jpm][1]));
            }
        }
        Eigenvalues::ApproxConstraint(eig_constraint) => {
            for jpm in 0..2 {
                let mut n = 0.0f64;
                for f in 0..F {
                    a[jpm] = a[jpm].max((fj[jpm + 1][f] - fj[jpm][f]).abs());
                    n += (vals[jpm + 2][f] - vals[jpm + 1][f]).powi(2);
                }
                if n == 0.0 {
                    a[jpm] = 0.0;
                } else {
                    a[jpm] /= n.sqrt();
                }
                a[jpm] = eig_constraint(t, tupm[jpm][0], a[jpm]).max(eig_constraint(
                    t,
                    tupm[jpm][1],
                    a[jpm],
                ));
            }
        }
    }
    let mut h: [[f64; F]; 2] = [[0.0; F]; 2];
    for jpm in 0..2 {
        for n in 0..F {
            h[jpm][n] =
                fpm[jpm][1][n] + fpm[jpm][0][n] - a[jpm] * (upm[jpm][1][n] - upm[jpm][0][n]);
        }
    }
    let mut res: [f64; F] = [0.0; F];
    for n in 0..F {
        res[n] = (h[1][n] - h[0][n]) / (2.0 * dx);
    }

    res
}