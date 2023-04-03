use std::{
    fs::File,
    io::{BufReader, Read},
};

use byteorder::{ByteOrder, LittleEndian};

use super::{HydroOutput, FREESTREAM_2D};

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

pub fn converge<const VX: usize, const VY: usize, const F: usize, const C: usize>(
    mut dt: f64,
    dtmin: f64,
    erpow: i32,
    fun: impl Fn(f64, f64) -> HydroOutput<VX, VY, F, C>,
) -> Option<()> {
    let mut ermin = dtmin.powi(erpow);
    let dtmul = 0.5;
    if ermin < 1e-15 {
        eprintln!(
            "dtmin^{}<1e-15 in converge, might not converge, set to 1e-15 for safety.",
            erpow
        );
        ermin = 1e-15;
    }
    let mut er = dt.powi(erpow);
    let mut f = fun(dt, er)?.0;
    println!("error convergence:");
    let update = |dt: f64| {
        let dt = dt * dtmul;
        let er = dt.powi(erpow);
        (dt, er)
    };
    (dt, er) = update(dt);

    while er > ermin {
        let f2 = fun(dt, er)?.0;
        let (ma, av) = compare(0, &f.0, &f2.0);
        println!("er: {:.3e}, max: {:.3e}, average: {:.3e}", er, ma, av);
        f = f2;
        (dt, er) = update(dt);
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

pub fn prepare_trento<const V: usize>(nb_trento: usize) -> Vec<[[f64; V]; V]> {
    let mut trentos = vec![[[0.0f64; V]; V]; nb_trento];
    let width = 1 + nb_trento.ilog10() as usize;
    for i in 0..nb_trento {
        trentos[i] = load_matrix(&format!("s{}/{:0>width$}.dat", V, i)).expect(&format!(
            "Could not load trento initial condition file \"e{V}/{i:0>width$}.dat\"."
        ));
    }
    trentos
}

pub fn prepare_trento_freestream<const V: usize>(
    nb_trento: usize,
) -> Vec<[[[f64; FREESTREAM_2D]; V]; V]> {
    let mut trentos = vec![[[[0.0f64; FREESTREAM_2D]; V]; V]; nb_trento];
    let width = 1 + nb_trento.ilog10() as usize;
    for tr in 0..nb_trento {
        let err = format!(
            "Could not load freestream trento initial condition file \"e{V}/{tr:0>width$}.dat\"."
        );
        let size = FREESTREAM_2D * std::mem::size_of::<f64>();
        let mut bytes = vec![0; size * V * V];
        let file = File::open(&format!("s{V}/data{tr:0>width$}.dat")).expect(&err);
        let mut buf = BufReader::new(file);
        buf.read_exact(&mut bytes).expect(&err);
        for y in 0..V {
            for x in 0..V {
                let i = x + y * V;
                LittleEndian::read_f64_into(
                    &bytes[i * size..(i + 1) * size],
                    &mut trentos[tr][y][x],
                );
            }
        }
    }
    trentos
}
