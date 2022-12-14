use sparse21::Matrix;
type Boundary<'a> = &'a dyn Fn(i32, usize) -> usize;
type Fun<'a> = &'a dyn Fn(&[f64], Boundary<'a>, i32) -> f64;

fn newton<'a>(
    fs: &'a Vec<(Fun<'a>, Boundary<'a>, i32)>,
    v: Vec<f64>,
    dt: f64,
    r: &Vec<Vec<f64>>,
) -> (Vec<f64>, usize) {
    let n = v.len();
    let nf = fs.len();
    assert!(
        n % nf == 0,
        "the size of the vector must be a multiple of the number of functions."
    );
    let nd = n / nf;
    let nr = r.len();
    let mut k = vec![0.0; nr * n];
    let ff = |vdtk: &Vec<f64>, k: &Vec<f64>, p: i32, all: bool| {
        let mut tmp = vec![0.0; n * nr];
        for i in 0..nr {
            for j in 0..nf {
                let ids: Vec<usize> = if all {
                    (0..nd).collect()
                } else {
                    (-fs[j].2..=fs[j].2).map(|l| fs[j].1(l + p, nd)).collect()
                };
                for l in ids {
                    let id = l + nd * (j + nf * i);
                    tmp[id] = fs[j].0(
                        &vdtk[(i * nf + j) * nd..(i * nf + j + 1) * nd],
                        fs[j].1,
                        l as i32,
                    ) - k[id];
                }
            }
        }
        tmp
    };
    let mut err = 1.0;
    let mut iterations = 0;
    while err > dt {
        let e = dt / nr as f64; // devide by 2 because 2 stages
        iterations += 1;
        let mut vdtk = (0..nr * n)
            .map(|i| {
                let j = i % n;
                let ni = i / n;
                v[j] + dt
                    * r[ni]
                        .iter()
                        .enumerate()
                        .fold(0.0, |acc, (ri, v)| acc + v * k[j + n * ri])
            })
            .collect::<Vec<_>>();
        let fu = ff(&vdtk, &k, 0, true);
        let mut m = Matrix::new();
        let mut count = 0;
        for i in 0..nr * n {
            let u = vdtk[i];
            let uk = k[i];
            let ni = i / n;
            let c = r[ni][ni];
            vdtk[i] = u + dt * c * e; // temporirilly change value for derivative computation
            k[i] = uk + e;
            let fpu = ff(&vdtk, &k, (i / (nr * nf)) as i32, false);
            for ir in 0..nr {
                for j in 0..nf {
                    for l in (-fs[j].2..=fs[j].2).map(|l| fs[j].1(l + (i / (nr * nf)) as i32, nd)) {
                        let id = l + nd * (j + nf * ir);
                        let d = (fpu[id] - fu[id]) / e;
                        if d != 0.0 {
                            m.add_element(id, i, d);
                            count += 1;
                        }
                    }
                }
            }
            vdtk[i] = u; // reset original value
            k[i] = uk;
        }
        if count > 0 {
            let mdks = m
                .solve(fu)
                .expect("Solving sparse system failed in Newton."); // the m stands for '-' as it is the opposite of the delta
            err = mdks.iter().fold(0.0f64, |acc, i| acc.max(i.abs()));
            k.iter_mut()
                .zip(mdks.iter())
                .for_each(|(k, mdk)| *k -= *mdk);
        } else {
            panic!("Nothing to do");
        }
    }
    (
        (0..n)
            .map(|i| {
                v[i] + dt
                    * r[r.len() - 1]
                        .iter()
                        .enumerate()
                        .fold(0.0, |acc, (ri, v)| acc + v * k[i + n * ri])
            })
            .collect::<Vec<_>>(),
        iterations,
    )
}

fn periodic(j: i32, n: usize) -> usize {
    if j >= n as i32 {
        j as usize - n
    } else if j < 0 {
        (n as i32 + j) as usize
    } else {
        j as usize
    }
}
fn flux(v: &[f64], boundary: Boundary, i: i32) -> f64 {
    let ivdx = 10.0f64;
    let n = v.len();
    ivdx.powf(2.0)
        * (-v[boundary(i - 2, n)] + 16.0 * v[boundary(i - 1, n)] - 30.0 * v[boundary(i, n)]
            + 16.0 * v[boundary(i + 1, n)]
            - v[boundary(i + 2, n)])
        / 12.0
}

fn minmod(theta: f64, a: f64, b: f64, c: f64) -> f64 {
    let minmod2 = |a: f64, b: f64| (a.signum() + b.signum()) / 2.0 * a.abs().min(b.abs());
    minmod2(theta * (a - b), minmod2((a - c) / 2.0, theta * (b - c)))
}

fn main() {
    let n = 100;
    let gaussian: Vec<f64> = (0..n)
        .map(|i| 100.0 * (-((i - n / 2) as f64).abs()).exp())
        .collect();
    let tot = gaussian.iter().fold(0.0, |acc, i| acc + i);
    let vals = [gaussian];
    let mut v: Vec<f64> = vals.into_iter().flat_map(|v| v.into_iter()).collect();
    let dt = 0.1;
    //let r = vec![vec![5.0 / 12.0, -1.0 / 12.0], vec![3.0 / 4.0, 1.0 / 4.0]];
    let r = vec![vec![1.0]];
    let t = 100;
    let f: Vec<(Fun, Boundary, i32)> = vec![(&flux, &periodic, 4)];
    let mul = f.iter().fold(0, |acc, (_, _, i)| acc.max(*i)) * 2 + 1;
    let tsteps = (t as f64 / dt) as usize;
    let mut tot_f = 0;
    for i in 0..tsteps {
        let (nv, c) = newton(&f, v, dt, &r);
        let c = c * (mul as usize + 1); // *(mul+1) because newton needs mul=2*'max distance interaction'+1 for the finite difference plus 1 for the reference
        v = nv;
        tot_f += c;
        if i % (1 + tsteps / 10) == 0 {
            println!("i: {}, c: {}", i, c);
        }
    }
    println!("v: {:?}", v);
    println!("expected mean: {}", tot / n as f64);
    println!(
        "tot_f: {}, tsteps: {}, ratio: {}",
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
}
