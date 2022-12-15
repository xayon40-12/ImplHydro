use crate::{newton::Boundary, utils::flux_limiter};

type Flux<'a, const C: usize> = &'a dyn Fn([f64; C]) -> f64;
type Constraints<'a, const F: usize, const C: usize> = &'a dyn Fn([f64; F]) -> [f64; C];
type Eigenvalues<'a, const C: usize> = &'a dyn Fn([f64; C]) -> f64;

pub enum Dir {
    X,
    Y,
}

pub fn kt<const F: usize, const VX: usize, const VY: usize, const C: usize, const N: usize>(
    vars: &[[[f64; F]; VX]; VY],
    [boundx, boundy]: &[Boundary; 2],
    [x, y]: [i32; 2],
    dir: Dir,
    flux: [Flux<C>; N],
    constraints: Constraints<F, C>,
    eigenvalues: Eigenvalues<C>,
    dx: f64,
    theta: f64,
) -> [f64; N] {
    let mut vals: [[f64; C]; 5] = [[0.0; C]; 5];
    for l in 0..5 {
        let (lx, ly) = match dir {
            Dir::X => (l as i32 - 2, 0),
            Dir::Y => (0, l as i32 - 2),
        };
        vals[l] = constraints(vars[boundy(y + ly, VY)][boundx(x + lx, VX)])
    }
    let mut deriv: [[f64; C]; 3] = [[0.0; C]; 3];
    for l in 0..3 {
        for c in 0..C {
            deriv[l][c] = flux_limiter(theta, vals[l][c], vals[l + 1][c], vals[l + 2][c]);
        }
    }
    let mut upm: [[[f64; C]; 2]; 2] = [[[0.0; C]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            let s = pm as f64 - 0.5;
            let j = jpm + pm;
            for c in 0..C {
                upm[jpm][pm][c] = vals[1 + j][c] + s * deriv[j][c];
            }
        }
    }
    let mut fpm: [[[f64; N]; 2]; 2] = [[[0.0; N]; 2]; 2];
    for jpm in 0..2 {
        for pm in 0..2 {
            for n in 0..N {
                fpm[jpm][pm][n] = flux[n](upm[jpm][pm]);
            }
        }
    }
    let mut a: [f64; 2] = [0.0; 2];
    for jpm in 0..2 {
        a[jpm] = eigenvalues(upm[jpm][0]).max(eigenvalues(upm[jpm][1]));
    }
    let mut h: [[f64; N]; 2] = [[0.0; N]; 2];
    for jpm in 0..2 {
        for n in 0..N {
            h[jpm][n] =
                fpm[jpm][1][n] + fpm[jpm][0][n] - a[jpm] * (upm[jpm][1][n] - upm[jpm][0][n]);
        }
    }
    let mut res: [f64; N] = [0.0; N];
    for n in 0..N {
        res[n] = (h[1][n] - h[0][n]) / (2.0 * dx);
    }

    res
}