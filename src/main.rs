use sparse21::Matrix;

fn flux(v: &[f64]) -> Vec<f64> {
    let ivdx = 10.0f64;
    let n = v.len();
    let periodic = |j: usize| {
        if j >= n {
            j - n
        } else {
            j
        }
    };
    return (0..n)
        .map(|i| {
            ivdx.powf(2.0)
                * (-v[periodic(i + n - 2)] + 16.0 * v[periodic(i + n - 1)] - 30.0 * v[i]
                    + 16.0 * v[periodic(i + 1)]
                    - v[periodic(i + 2)])
                / 12.0
        })
        .collect();
}

fn exp(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    (0..n).map(|i| -v[i]).collect()
}

fn newton(
    f: &dyn Fn(&[f64]) -> Vec<f64>,
    v: Vec<f64>,
    dt: f64,
    r: &Vec<Vec<f64>>,
) -> (Vec<f64>, usize) {
    let n = v.len();
    let rn = r.len();
    let mut k = vec![0.0; rn * n];
    let ff = |vdtk: &Vec<f64>, k: &Vec<f64>| {
        (0..rn)
            .fold(vec![].into_iter(), |acc, i| {
                acc.chain(f(&vdtk[i * n..(i + 1) * n]).into_iter())
                    .collect::<Vec<_>>()
                    .into_iter()
            })
            .zip(k.iter())
            .map(|(v, k)| v - k)
            .collect::<Vec<_>>()
    };
    let mut err = 1.0;
    let mut iterations = 0;
    while err > dt {
        let e = dt / rn as f64; // devide by 2 because 2 stages
        iterations += 1;
        let mut vdtk = (0..rn * n)
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
        let fu = ff(&vdtk, &k);
        let mut m = Matrix::new();
        let mut count = 0;
        for i in 0..rn * n {
            let u = vdtk[i];
            let uk = k[i];
            let ni = i / n;
            let c = r[ni][ni];
            vdtk[i] = u + dt * c * e; // temporirilly change value for derivative computation
            k[i] = uk + e;
            for (j, fpu) in ff(&vdtk, &k).iter().enumerate() {
                let d = (fpu - fu[j]) / e;
                if d != 0.0 {
                    //println!("i: {}, j: {}, d: {}", i, j, d);
                    m.add_element(j, i, d);
                    count += 1;
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

fn main() {
    let n = 100;
    let mut v: Vec<f64> = (0..n)
        .map(|i| 100.0 * (-((i - n / 2) as f64).abs()).exp())
        .collect();
    let dt = 0.1;
    //let r = vec![vec![5.0 / 12.0, -1.0 / 12.0], vec![3.0 / 4.0, 1.0 / 4.0]];
    let r = vec![vec![1.0]];
    let t = 100;
    let f = &flux;
    let tsteps = (t as f64 / dt) as usize;
    let mut tot_iter = 0;
    for i in 0..tsteps {
        let (nv, c) = newton(f, v, dt, &r);
        v = nv;
        tot_iter += c;
        if i % (1 + tsteps / 10) == 0 {
            println!("i: {}, c: {}", i, c);
        }
    }
    let tot_f = tot_iter * 6; // 1 f + 3 f for space derivatives
                              //println!("exp(-1): {}", (-1.0f64).exp());
    println!("v: {:?}", v);
    println!(
        "tot_iter: {}, tot_f: {}, tsteps: {}, ratio: {}",
        tot_iter,
        tot_f,
        tsteps,
        tot_f as f64 / tsteps as f64
    );
}
