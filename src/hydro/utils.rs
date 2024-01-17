use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read},
};

use crate::solver::context::{Arr, BArr, DIM};
use boxarray::boxarray;
use byteorder::{ByteOrder, LittleEndian};

use super::{HydroOutput, FREESTREAM_2D};

#[derive(Clone)]
pub enum Coordinate {
    Cartesian,
    Milne,
}

pub fn eigenvaluesk(dpde: f64, ut: f64, uk: f64) -> f64 {
    let vs2 = dpde;
    let a = ut * uk * (1.0 - vs2);
    let b = (ut * ut - uk * uk - (ut * ut - uk * uk - 1.0) * vs2) * vs2;
    let d = ut * ut - (ut * ut - 1.0) * vs2;
    if b < 0.0 || d == 0.0 {
        panic!("\nvs2:\na: {a}, b:{b}, d:{d}\n");
    }
    (a.abs() + b.sqrt()) / d
}

pub fn compare<const VX: usize, const VY: usize, const VZ: usize, const F: usize>(
    f: usize,
    v1: &Arr<F, VX, VY, VZ>,
    v2: &Arr<F, VX, VY, VZ>,
) -> (f64, f64) {
    let mut average: f64 = 0.0;
    let mut maximum: f64 = 0.0;
    for vz in 0..VZ {
        for vy in 0..VY {
            for vx in 0..VX {
                let v1 = v1[vz][vy][vx][f];
                let v2 = v2[vz][vy][vx][f];
                let d = (v1 - v2).abs() / v1.abs().max(v2.abs());
                average += d;
                maximum = maximum.max(d);
            }
        }
    }
    average /= (VX * VY * VZ) as f64;

    (maximum, average)
}

pub fn converge<
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const F: usize,
    const C: usize,
>(
    mut dt: f64,
    dtmin: f64,
    fun: impl Fn(f64) -> HydroOutput<VX, VY, VZ, F, C>,
) -> Option<()> {
    let dtmul = 0.5;
    let update = |dt: f64| {
        let dt = dt * dtmul;
        dt
    };
    let mut f: Option<(BArr<F, VX, VY, VZ>, BArr<C, VX, VY, VZ>)> = None;
    // let mut f = fun(dt)?.0;
    // println!("error convergence:");
    // dt = update(dt);

    while dt > dtmin {
        if let Some(f2) = fun(dt) {
            let f2 = f2.0;
            if let Some(f) = f {
                let (ma, av) = compare(0, &f.0, &f2.0);
                println!("dt: {:.3e}, max: {:.3e}, average: {:.3e}", dt, ma, av);
            }
            f = Some(f2);
        }
        dt = update(dt);
    }
    println!("");

    Some(())
}

pub fn load_matrix_2d<const VX: usize, const VY: usize>(
    filename: &str,
) -> std::io::Result<Box<[[f64; VX]; VY]>> {
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
            "Wrong matrix size in load_matrix_3d: expected {}x{}, found {}x{}",
            VX, VY, vx, vy
        );
    }
    let mut mat: Box<[[f64; VX]; VY]> = boxarray(0.0);
    for j in 0..VY {
        for i in 0..VX {
            mat[j][i] = arr[j][i];
        }
    }

    Ok(mat)
}

pub fn load_matrix_3d<const VX: usize, const VY: usize, const VZ: usize>(
    filename: &str,
) -> std::io::Result<Box<[[[f64; VX]; VY]; VZ]>> {
    let mut file = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let arr = contents
        .split("\n")
        .filter(|l| !l.starts_with("#") && l.len() > 0)
        .map(|l| {
            l.trim()
                .split(" ")
                .map(|v| {
                    v.parse::<f64>().expect(&format!(
                        "Cannot parse f64 in load_matrix call for file \"{}\"",
                        filename
                    ))
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<Vec<_>>>();

    let vyz = arr.len();
    let vz = vyz / VY;
    let vy = VY;
    let vx = arr[0].len();

    if vz % 2 != VZ % 2 || VZ > vz {
        panic!("Wrong matrix size in load_matrix_3d: the number of cells in z/eta direction must be smaller or equal to the loaded data and they should both be the same even or odd, expected {} found {}.", vz, VZ);
    }
    if vx != VX || vy != VY {
        panic!(
            "Wrong matrix size in load_matrix_3d: expected {}x{}x{}, found {}x{}x{} (the third component is allowed to be smaller)",
            VX, VY, VZ, vx, vy, vz
        );
    }
    let mut mat: Box<[[[f64; VX]; VY]; VZ]> = boxarray(0.0);
    let sk = (vz - VZ) / 2;
    for k in 0..VZ {
        for j in 0..VY {
            for i in 0..VX {
                mat[k][j][i] = arr[j + VY * (sk + k)][i];
            }
        }
    }

    Ok(mat)
}

pub fn prepare_trento_2d<const V: usize>(nb_trento: usize) -> Vec<Box<[[f64; V]; V]>> {
    let mut trentos: Vec<Box<[[f64; V]; V]>> = vec![boxarray(0.0); nb_trento];
    let width = 1 + (nb_trento - 1).max(1).ilog10() as usize;
    for i in 0..nb_trento {
        trentos[i] = load_matrix_2d(&format!("s{}/{:0>width$}.dat", V, i)).expect(&format!(
            "Could not load trento initial condition file \"s{V}/{i:0>width$}.dat\"."
        ));
    }
    trentos
}

pub fn prepare_trento_3d<const XY: usize, const Z: usize>(
    nb_trento: usize,
) -> (Vec<Box<[[[f64; XY]; XY]; Z]>>, Option<[f64; DIM]>) {
    let mut trentos: Vec<Box<[[[f64; XY]; XY]; Z]>> = vec![boxarray(0.0); nb_trento];
    let mut dxs = None;
    let width = 1 + (nb_trento - 1).max(1).ilog10() as usize;
    for i in 0..nb_trento {
        const LOAD2D: bool = false;
        // const LOAD2D: bool = true;
        if LOAD2D {
            let trento_2d =
                load_matrix_2d(&format!("s{}/{:0>width$}.dat", XY, i)).expect(&format!(
                    "Could not load trento initial condition file \"s{XY}/{i:0>width$}.dat\"."
                ));
            for z in 0..Z {
                trentos[i][z] = *trento_2d;
            }
        } else {
            trentos[i] = load_matrix_3d(&format!("s{}/{:0>width$}.dat", XY, i)).expect(&format!(
                "Could not load trento initial condition file \"s{XY}/{i:0>width$}.dat\"."
            ));
            if let Ok(str) = std::fs::read_to_string(&format!("s{XY}/info.txt")) {
                let m: HashMap<String, f64> = str
                    .trim()
                    .split("\n")
                    .map(|s| {
                        let v = s.split(": ").collect::<Vec<_>>();
                        (
                            v[0].to_string(),
                            v[1].parse()
                                .expect("Could not parse lattice spacing for trento info.txt."),
                        )
                    })
                    .collect();
                dxs = Some([m["dx"], m["dy"], m["detas"]]);
            }
        }
    }
    (trentos, dxs)
}

pub fn prepare_trento_2d_freestream<const V: usize>(
    nb_trento: usize,
) -> Vec<Arr<FREESTREAM_2D, V, V, 1>> {
    let mut trentos = vec![[[[[0.0f64; FREESTREAM_2D]; V]; V]; 1]; nb_trento];
    let width = 1 + (nb_trento - 1).max(1).ilog10() as usize;
    for tr in 0..nb_trento {
        let err = format!(
            "Could not load freestream trento initial condition file \"s{V}/{tr:0>width$}.dat\"."
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
                    &mut trentos[tr][0][y][x],
                );
            }
        }
    }
    trentos
}
