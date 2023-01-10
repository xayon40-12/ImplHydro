use std::{fs::File, io::Read};

use super::{hydro1d, hydro2d, Pressure};

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

pub fn init_from_energy_1d<'a, const VX: usize>(
    t0: f64,
    es: [f64; VX],
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn(usize, f64) -> [f64; 2] + 'a> {
    Box::new(move |i, _| {
        // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
        let e = es[i].max(1e-100);
        let vars = [e, p(e), dpde(e), 1.0, 0.0];
        [hydro1d::f00(t0, vars), hydro1d::f01(t0, vars)]
    })
}

pub fn init_from_energy_2d<'a, const VX: usize, const VY: usize>(
    t0: f64,
    es: [[f64; VX]; VY],
    p: Pressure<'a>,
    dpde: Pressure<'a>,
) -> Box<dyn Fn((usize, usize), (f64, f64)) -> [f64; 3] + 'a> {
    Box::new(move |(i, j), _| {
        // let e = if x == 0.0 && y == 0.0 { 10.0 } else { 1e-100 };
        let e = es[j][i].max(1e-100);
        let vars = [e, p(e), dpde(e), 1.0, 0.0, 0.0];
        [
            hydro2d::f00(t0, vars),
            hydro2d::f01(t0, vars),
            hydro2d::f02(t0, vars),
        ]
    })
}
