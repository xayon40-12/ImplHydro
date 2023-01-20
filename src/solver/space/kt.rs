use crate::solver::{context::Boundary, utils::flux_limiter, Transform};

use super::{Dir, Eigenvalues, Flux};

pub fn kt<const F: usize, const VX: usize, const VY: usize, const C: usize>(
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
    theta: f64,
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
    let mut deriv: [[f64; F]; 3] = [[0.0; F]; 3];
    let mut fj: [[f64; F]; 3] = [[0.0; F]; 3];
    for l in 0..3 {
        for f in 0..F {
            deriv[l][f] = flux_limiter(theta, vals[l][f], vals[l + 1][f], vals[l + 2][f]);
            fj[l][f] = flux[f](t, tvals[l + 1]);
        }
    }
    let mut upm: [[[f64; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
    let mut tupm: [[[f64; C]; 2]; 2] = [[[0.0; C]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            let s = pm as f64 - 0.5;
            let j = jpm + pm;
            for f in 0..F {
                upm[jpm][pm][f] = vals[1 + j][f] - s * deriv[j][f];
            }
            upm[jpm][pm] = constraints(t, post_flux_limiter(t, upm[jpm][pm]));
            tupm[jpm][pm] = transform(t, upm[jpm][pm]);
        }
    }
    let mut fpm: [[[f64; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
    for n in 0..F {
        for jpm in 0..2 {
            for pm in 0..2 {
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
