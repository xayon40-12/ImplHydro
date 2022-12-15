use itertools::Itertools;
use sparse21::Matrix;

pub type Boundary<'a> = &'a dyn Fn(i32, usize) -> usize;
pub type Fun<'a, const F: usize, const VX: usize, const VY: usize> =
    &'a dyn Fn(&[[[f64; F]; VX]; VY], &[Boundary; 2], [i32; 2]) -> [f64; F];

pub struct Context<'a, 'b, const F: usize, const VX: usize, const VY: usize, const S: usize> {
    pub fun: Fun<'a, F, VX, VY>,
    pub boundary: &'b [Boundary<'b>; 2],
    pub local_interaction: [i32; 2],
    pub vs: [[[f64; VX]; VY]; F],
    pub k: [[[[f64; F]; VX]; VY]; S],
    pub integrated: [bool; F],
    pub r: [[f64; S]; S],
    pub dt: f64,
    pub er: f64,
    pub tbeg: f64,
    pub tend: f64,
}

fn sub<const F: usize>(a: [f64; F], b: [f64; F]) -> [f64; F] {
    let mut res = [0.0; F];
    for i in 0..F {
        res[i] = a[i] - b[i];
    }
    res
}

fn flat<const F: usize, const VX: usize, const VY: usize, const S: usize>(
    a: [[[[f64; F]; VX]; VY]; S],
) -> Vec<f64> {
    let mut res = vec![0.0; F * VX * VY * S];
    for s in 0..S {
        for f in 0..F {
            for vy in 0..VY {
                for vx in 0..VX {
                    let x = f + F * (vx + VX * (vy + VY * s));
                    res[x] = a[s][vy][vx][f];
                }
            }
        }
    }
    res
}

pub fn newton<const F: usize, const VX: usize, const VY: usize, const S: usize>(
    Context {
        fun,
        boundary,
        local_interaction,
        vs,
        k,
        integrated,
        r,
        dt,
        er,
        tbeg: _,
        tend: _,
    }: &mut Context<F, VX, VY, S>,
) -> usize {
    let [sizex, sizey] = *local_interaction;
    let ff = |vdtk: &[[[[f64; F]; VX]; VY]; S],
              k: &[[[[f64; F]; VX]; VY]; S],
              [px, py]: [i32; 2],
              all: bool| {
        let mut tmp = [[[[0.0; F]; VX]; VY]; S];
        for s in 0..S {
            let (idsx, idsy): (Vec<usize>, Vec<usize>) = if all {
                ((0..VX).collect(), (0..VY).collect())
            } else {
                (
                    (-sizex..=sizex)
                        .map(|l| boundary[0](l + px, VX))
                        .unique()
                        .collect(),
                    (-sizey..=sizey)
                        .map(|l| boundary[1](l + py, VY))
                        .unique()
                        .collect(),
                )
            };
            for vy in idsy {
                for vx in &idsx {
                    tmp[s][vy][*vx] = sub(
                        fun(&vdtk[s], boundary, [*vx as i32, vy as i32]),
                        k[s][vy][*vx],
                    );
                }
            }
        }
        tmp
    };
    let mut err = 1.0;
    let mut iterations = 0;
    while err > *er {
        let e = *dt / S as f64; // devide by S because S stages
        iterations += 1;
        let mut vdtk = [[[[0.0; F]; VX]; VY]; S];
        for f in 0..F {
            for vy in 0..VY {
                for vx in 0..VX {
                    for s in 0..S {
                        if integrated[f] {
                            for s1 in 0..S {
                                vdtk[s][vy][vx][f] +=
                                    vs[f][vy][vx] + *dt * r[s][s1] * k[s1][vy][vx][f];
                            }
                        } else {
                            vdtk[s][vy][vx][f] = k[s][vy][vx][f];
                        }
                    }
                }
            }
        }

        let fu = ff(&vdtk, &k, [0, 0], true);
        let mut m = Matrix::new();
        let mut count = 0;
        for s0 in 0..S {
            for f0 in 0..F {
                for vy0 in 0..VY {
                    for vx0 in 0..VX {
                        let u = vdtk[s0][vy0][vx0][f0];
                        let uk = k[s0][vy0][vx0][f0];
                        let c = r[s0][s0];
                        if integrated[f0] {
                            vdtk[s0][vy0][vx0][f0] = u + *dt * c * e;
                        } else {
                            vdtk[s0][vy0][vx0][f0] = u + e;
                        }
                        k[s0][vy0][vx0][f0] = uk + e;

                        let fpu = ff(&vdtk, &k, [vx0 as i32, vy0 as i32], false);
                        for s1 in 0..S {
                            for f1 in 0..F {
                                for vy1 in (-sizey..=sizey)
                                    .map(|l| boundary[1](l + vy0 as i32, VY))
                                    .unique()
                                {
                                    for vx1 in (-sizex..=sizex)
                                        .map(|l| boundary[0](l + vx0 as i32, VX))
                                        .unique()
                                    {
                                        let d = (fpu[s1][vy1][vx1][f1] - fu[s1][vy1][vx1][f1]) / e;
                                        if d != 0.0 {
                                            let x = f0 + F * (vy0 + VY * (vx0 + VX * s0));
                                            let y = f1 + F * (vy1 + VY * (vx1 + VX * s1));
                                            m.add_element(y, x, d);
                                            count += 1;
                                        }
                                    }
                                }
                            }
                        }

                        vdtk[s0][vy0][vx0][f0] = u;
                        k[s0][vy0][vx0][f0] = uk;
                    }
                }
            }
        }
        if count > 0 {
            let mdks = m
                .solve(flat(fu))
                .expect("Solving sparse system failed in Newton."); // the m stands for '-' as it is the opposite of the delta
            err = mdks.iter().fold(0.0f64, |acc, i| acc.max(i.abs()));
            for f in 0..F {
                for vy in 0..VY {
                    for vx in 0..VX {
                        for s in 0..S {
                            let x = f + F * (vy + VY * (vx + VX * s));
                            k[s][vy][vx][f] -= mdks[x];
                        }
                    }
                }
            }
        } else {
            panic!("Nothing to do");
        }
    }
    for f in 0..F {
        for vy in 0..VY {
            for vx in 0..VX {
                if integrated[f] {
                    for s in 0..S {
                        vs[f][vy][vx] += *dt * r[S - 1][s] * k[s][vy][vx][f];
                    }
                } else {
                    vs[f][vy][vx] = k[S - 1][vy][vx][f];
                }
            }
        }
    }
    iterations
}
