use std::{fs::File, io::Read};

use super::HydroOutput;

pub fn compare<const VX: usize, const VY: usize, const F: usize>(
    f: usize,
    v1: &[[[f64; F]; VX]; VY],
    v2: &[[[f64; F]; VX]; VY],
) -> (f64, f64) {
    let mut average: f64 = 0.0;
    let mut maximum: f64 = 0.0;
    for j in 0..VY {
        for i in 0..VX {
            let v1 = v1[j][i][f];
            let v2 = v2[j][i][f];
            let d = (v1 - v2).abs() / v1.abs().max(v2.abs());
            average += d;
            maximum = maximum.max(d);
        }
    }
    average /= (VX * VY) as f64;

    (maximum, average)
}

pub fn converge<const VX: usize, const VY: usize, const F: usize>(
    mut er: f64,
    mut ermin: f64,
    fun: impl Fn(f64) -> HydroOutput<VX, VY, F>,
) -> Option<()> {
    if ermin < 1e-15 {
        eprintln!("ermin<1e-15 in converge, might not converge, set to ermin=1e-15 for safety.");
        ermin = 1e-15;
    }
    let mut f = fun(er)?.0;
    println!("error convergence:");
    er *= 0.1;
    while er > ermin {
        let f2 = fun(er)?.0;
        let (ma, av) = compare(0, &f, &f2);
        println!("er: {:.3e}, max: {:.3e}, average: {:.3e}", er, ma, av);
        f = f2;
        er *= 0.1;
    }
    println!("");

    Some(())
}

pub fn load_matrix<const VX: usize, const VY: usize>(
    filename: &str,
) -> std::io::Result<[[f64; VX]; VY]> {
    let mut file = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let arr = contents
        .split("\n")
        .filter(|l| !l.starts_with("#") && l.len() > 0)
        .map(|l| {
            l.split(" ")
                .map(|v| {
                    v.parse::<f64>().expect(&format!(
                        "Cannot parse f64 in load_matrix call for file \"{}\"",
                        filename
                    ))
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<Vec<_>>>();

    let vx = arr[0].len();
    let vy = arr.len();
    if vy != VY || vx != VX {
        panic!(
            "Wrong matrix size in load_matrix: expected {}x{}, found {}x{}",
            VX, VY, vx, vy
        );
    }
    let mut mat = [[0.0f64; VX]; VY];
    for j in 0..VY {
        for i in 0..VX {
            mat[j][i] = arr[j][i];
        }
    }

    Ok(mat)
}

pub fn prepare_trento<const V: usize, const TRENTO: usize>() -> [[[f64; V]; V]; TRENTO] {
    let mut trentos = [[[0.0f64; V]; V]; TRENTO];
    for i in 0..TRENTO {
        trentos[i] = load_matrix(&format!("e{}/{:0>2}.dat", V, i)).expect(&format!(
            "Could not load trento initial condition file \"e{}/{:0>2}.dat\".",
            V, i
        ));
    }
    trentos
}
