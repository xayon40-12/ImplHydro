use sparse21::Matrix;
pub type Boundary<'a> = &'a dyn Fn(i32, usize) -> usize;
pub type Fun<'a, const F: usize, const V: usize> =
    &'a dyn Fn(&[[f64; F]; V], Boundary<'a>, i32) -> [f64; F];

fn sub<const F: usize>(a: [f64; F], b: [f64; F]) -> [f64; F] {
    let mut res = [0.0; F];
    for i in 0..F {
        res[i] = a[i] - b[i];
    }
    res
}

fn flat<const F: usize, const V: usize, const S: usize>(a: [[[f64; F]; V]; S]) -> Vec<f64> {
    let mut res = vec![0.0; F * V * S];
    for s in 0..S {
        for f in 0..F {
            for v in 0..V {
                let x = f + F * (v + V * s);
                res[x] += a[s][v][f];
            }
        }
    }
    res
}

pub fn newton<'a, const F: usize, const V: usize, const S: usize>(
    (fun, boundary, size): (Fun<'a, F, V>, Boundary<'a>, i32),
    mut vs: [[f64; V]; F],
    dt: f64,
    r: &[[f64; S]; S],
) -> ([[f64; V]; F], usize) {
    let mut k = [[[0.0; F]; V]; S];
    let ff = |vdtk: &[[[f64; F]; V]; S], k: &[[[f64; F]; V]; S], p: i32, all: bool| {
        let mut tmp = [[[0.0; F]; V]; S];
        for s in 0..S {
            let ids: Vec<usize> = if all {
                (0..V).collect()
            } else {
                (-size..=size).map(|l| boundary(l + p, V)).collect()
            };
            for v in ids {
                tmp[s][v] = sub(fun(&vdtk[s], boundary, v as i32), k[s][v]);
            }
        }
        tmp
    };
    let mut err = 1.0;
    let mut iterations = 0;
    while err > dt {
        let e = dt / S as f64; // devide by 2 because 2 stages
        iterations += 1;
        let mut vdtk = [[[0.0; F]; V]; S];
        for f in 0..F {
            for v in 0..V {
                for s in 0..S {
                    for s1 in 0..S {
                        vdtk[s][v][f] += vs[f][v] + dt * r[s][s1] * k[s1][v][f];
                    }
                }
            }
        }

        let fu = ff(&vdtk, &k, 0, true);
        let mut m = Matrix::new();
        let mut count = 0;
        for s0 in 0..S {
            for f0 in 0..F {
                for v0 in 0..V {
                    let u = vdtk[s0][v0][f0];
                    let uk = k[s0][v0][f0];
                    let c = r[s0][s0];
                    vdtk[s0][v0][f0] = u + dt * c * e;
                    k[s0][v0][f0] = uk + e;

                    let fpu = ff(&vdtk, &k, v0 as i32, false);
                    for s1 in 0..S {
                        for f1 in 0..F {
                            for v1 in (-size..=size).map(|l| boundary(l + v0 as i32, V)) {
                                let d = (fpu[s1][v1][f1] - fu[s1][v1][f1]) / e;
                                if d != 0.0 {
                                    let x = f0 + F * (v0 + V * s0);
                                    let y = f1 + F * (v1 + V * s1);
                                    m.add_element(y, x, d);
                                    count += 1;
                                }
                            }
                        }
                    }

                    vdtk[s0][v0][f0] = u;
                    k[s0][v0][f0] = uk;
                }
            }
        }
        if count > 0 {
            let mdks = m
                .solve(flat(fu))
                .expect("Solving sparse system failed in Newton."); // the m stands for '-' as it is the opposite of the delta
            err = mdks.iter().fold(0.0f64, |acc, i| acc.max(i.abs()));
            for f in 0..F {
                for v in 0..V {
                    for s in 0..S {
                        let x = f + F * (v + V * s);
                        k[s][v][f] -= mdks[x];
                    }
                }
            }
        } else {
            panic!("Nothing to do");
        }
    }
    for f in 0..F {
        for v in 0..V {
            for s in 0..S {
                vs[f][v] += dt * r[S - 1][s] * k[s][v][f];
            }
        }
    }
    (vs, iterations)
}
