use crate::solver::context::{Arr, Boundary, DIM};
use crate::solver::{utils::flux_limiter, Constraint, Transform};
use crate::FLOAT;

use super::{Dir, Eigenvalues, FluxInfo, InDir};

pub fn kt<
    const N: usize,
    const F: usize,
    const C: usize,
    const SD: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
>(
    (vs, _trs): (&Arr<F, VX, VY, VZ>, &Arr<C, VX, VY, VZ>),
    bound: Boundary<F, VX, VY, VZ>,
    [x, y, z]: [i32; DIM],
    t: FLOAT,
    flux_infos: [InDir<FluxInfo<F, C, SD>>; N],
    constraints: Constraint<F, C>,
    pre_flux_limiter: Transform<F>,
    post_flux_limiter: Transform<F>,
    dxs: [FLOAT; DIM],
    theta: FLOAT,
) -> [([FLOAT; F], [FLOAT; SD]); N] {
    let mut res = [([0.0 as FLOAT; F], [0.0 as FLOAT; SD]); N];
    // let res: [InDir<([FLOAT; F], [FLOAT; SD])>; N];
    // for i in 0..flux_infos.len() {
    for (i, flux_info) in flux_infos.into_iter().enumerate() {
        let dir = flux_info.dir();
        let dx = match dir {
            Dir::X => dxs[0],
            Dir::Y => dxs[1],
            Dir::Z => dxs[2],
        };
        let FluxInfo {
            flux,
            secondary,
            eigenvalues,
        } = flux_info.iner();
        res[i] = {
            let mut vals: [[FLOAT; F]; 5] = [[0.0; F]; 5];
            let mut sd: [[FLOAT; SD]; 5] = [[0.0; SD]; 5];
            let mut pvals: [[FLOAT; F]; 5] = [[0.0; F]; 5];
            let mut tvals: [[FLOAT; C]; 5] = [[0.0; C]; 5];
            for l in 0..5 {
                let (lx, ly, lz) = match dir {
                    Dir::X => (l as i32 - 2, 0, 0),
                    Dir::Y => (0, l as i32 - 2, 0),
                    Dir::Z => (0, 0, l as i32 - 2),
                };
                (vals[l], tvals[l]) = constraints(t, bound([x + lx, y + ly, z + lz], &vs));
                sd[l] = secondary(t, tvals[l]);
                pvals[l] = pre_flux_limiter(t, vals[l]);
            }
            let mut deriv: [[FLOAT; F]; 3] = [[0.0; F]; 3];
            let mut fj: [[FLOAT; F]; 3] = [[0.0; F]; 3];
            for l in 0..3 {
                for f in 0..F {
                    deriv[l][f] =
                        flux_limiter(theta, pvals[l][f], pvals[l + 1][f], pvals[l + 2][f]);
                }
                fj[l] = flux(t, tvals[l + 1]);
            }
            let mut upm: [[[FLOAT; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
            let mut tupm: [[[FLOAT; C]; 2]; 2] = [[[0.0; C]; 2]; 2];
            for jpm in 0..2 {
                for pm in 0..2 {
                    let s = pm as FLOAT - 0.5;
                    let j = jpm + pm;
                    for f in 0..F {
                        upm[jpm][pm][f] = pvals[1 + j][f] - s * deriv[j][f];
                    }
                    upm[jpm][pm] = post_flux_limiter(t, upm[jpm][pm]);
                    (upm[jpm][pm], tupm[jpm][pm]) = constraints(t, upm[jpm][pm]);
                }
            }
            let mut fpm: [[[FLOAT; F]; 2]; 2] = [[[0.0; F]; 2]; 2];
            for jpm in 0..2 {
                for pm in 0..2 {
                    fpm[jpm][pm] = flux(t, tupm[jpm][pm]);
                }
            }
            let mut a: [FLOAT; 2] = [0.0; 2];
            match eigenvalues {
                Eigenvalues::Analytical(eig_analytical) => {
                    for jpm in 0..2 {
                        a[jpm] =
                            eig_analytical(t, tupm[jpm][0]).max(eig_analytical(t, tupm[jpm][1]));
                    }
                }
                Eigenvalues::ApproxConstraint(eig_constraint) => {
                    for jpm in 0..2 {
                        let mut n = 0.0 as FLOAT;
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
            let mut h: [[FLOAT; F]; 2] = [[0.0; F]; 2];
            let mut sdh: [[FLOAT; F]; 2] = [[0.0; F]; 2];
            let cubic = [-0.0625, 0.5625, 0.5625, -0.0625]; // reconstruct secondary at the boundary using cubic polynomial
            for jpm in 0..2 {
                for n in 0..F {
                    h[jpm][n] = fpm[jpm][1][n] + fpm[jpm][0][n]
                        - a[jpm] * (upm[jpm][1][n] - upm[jpm][0][n]);
                }
                for n in 0..SD {
                    for i in 0..4 {
                        sdh[jpm][n] += cubic[i] * sd[i + jpm][n];
                    }
                }
            }
            let mut res: [FLOAT; F] = [0.0; F];
            let mut sdres: [FLOAT; SD] = [0.0; SD];
            for n in 0..F {
                res[n] = (h[1][n] - h[0][n]) / (2.0 * dx);
            }
            for n in 0..SD {
                sdres[n] = (sdh[1][n] - sdh[0][n]) / dx;
            }

            (res, sdres)
        };
    }
    res
}
