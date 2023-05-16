use crate::solver::context::Boundary;
use crate::solver::{utils::flux_limiter, Constraint, Transform};

use super::{Dir, Eigenvalues, Flux};

pub fn kt<const F: usize, const VX: usize, const VY: usize, const C: usize, const SD: usize>(
    (vs, _trs): (&[[[f64; F]; VX]; VY], &[[[f64; C]; VX]; VY]),
    bound: Boundary<F, VX, VY>,
    [x, y]: [i32; 2],
    dir: Dir,
    t: f64,
    flux: Flux<F, C>,
    secondary: Flux<SD, C>,
    constraints: Constraint<F, C>,
    eigenvalues: Eigenvalues<C>,
    pre_flux_limiter: Transform<F>,
    post_flux_limiter: Transform<F>,
    dx: f64,
    theta: f64,
) -> ([f64; F], [f64; SD]) {
    let mut vals: [[f64; F]; 5] = [[0.0; F]; 5];
    let mut tvals: [[f64; C]; 5] = [[0.0; C]; 5];
    for l in 0..5 {
        let (lx, ly) = match dir {
            Dir::X => (l as i32 - 2, 0),
            Dir::Y => (0, l as i32 - 2),
        };
        (vals[l], tvals[l]) = constraints(t, bound([x + lx, y + ly], &vs));
        vals[l] = pre_flux_limiter(t, vals[l]);
    }
    let mut deriv: [[f64; F]; 3] = [[0.0; F]; 3];
    // let mut tderiv: [[f64; C]; 3] = [[0.0; C]; 3];
    let mut fj: [[f64; F]; 3] = [[0.0; F]; 3];
    for l in 0..3 {
        for f in 0..F {
            deriv[l][f] = flux_limiter(theta, vals[l][f], vals[l + 1][f], vals[l + 2][f]);
        }
        // for c in 0..C {
        //     tderiv[l][c] = flux_limiter(theta, tvals[l][c], tvals[l + 1][c], tvals[l + 2][c]);
        // }
        fj[l] = flux(t, tvals[l + 1]);
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
            // for c in 0..C {
            //     tupm[jpm][pm][c] = tvals[1 + j][c] - s * tderiv[j][c];
            // }
            upm[jpm][pm] = post_flux_limiter(t, upm[jpm][pm]);
            (upm[jpm][pm], tupm[jpm][pm]) = constraints(t, upm[jpm][pm]);
        }
    }
    let mut fpm: [[[f64; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
    let mut sdpm: [[[f64; SD]; 2]; 2] = [[[0.0; SD]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            fpm[jpm][pm] = flux(t, tupm[jpm][pm]);
            sdpm[jpm][pm] = secondary(t, tupm[jpm][pm]);
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
    let mut sdh: [[f64; F]; 2] = [[0.0; F]; 2];
    for jpm in 0..2 {
        for n in 0..F {
            h[jpm][n] =
                fpm[jpm][1][n] + fpm[jpm][0][n] - a[jpm] * (upm[jpm][1][n] - upm[jpm][0][n]);
        }
        for n in 0..SD {
            sdh[jpm][n] = sdpm[jpm][1][n] + sdpm[jpm][0][n];
        }
    }
    let mut res: [f64; F] = [0.0; F];
    let mut sdres: [f64; SD] = [0.0; SD];
    for n in 0..F {
        res[n] = (h[1][n] - h[0][n]) / (2.0 * dx);
    }
    for n in 0..SD {
        sdres[n] = (sdh[1][n] - sdh[0][n]) / (2.0 * dx);
    }

    // let m = res
    //     .iter()
    //     .map(|v| v.abs())
    //     .max_by(|a, b| a.total_cmp(b))
    //     .unwrap();
    // let msd = sdres
    //     .iter()
    //     .map(|v| v.abs())
    //     .max_by(|a, b| a.total_cmp(b))
    //     .unwrap();
    // let mm = 160.0;
    // if m > mm || msd > mm {
    //     println!("{:e} {:e}", m, msd);
    // }

    (res, sdres)
}
