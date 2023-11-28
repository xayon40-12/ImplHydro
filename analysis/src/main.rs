use num::Complex;
use std::{collections::HashMap, iter::repeat};

#[derive(Debug, Clone, Copy)]
pub struct Particle {
    pub eta: f64,
    pub pt: f64,
    pub phi: f64,
}

#[derive(Debug, Clone)]
pub struct FreezeOut {
    pub particles: Vec<Particle>,
}

#[derive(Debug, Clone)]
pub struct Event {
    pub freezeouts: Vec<FreezeOut>,
}

impl Event {
    pub fn multiplicity(&self) -> usize {
        self.freezeouts
            .iter()
            .map(|f| f.particles.len())
            .reduce(|a, b| a + b)
            .unwrap()
            / self.freezeouts.len()
    }
}

// def jackknife(fn, ns, ds):
//    n = sum(ns)
//    d = sum(ds)
//    value = fn(n, d)
//    l = len(ns)
//    if l > 1:
//       fs = [fn(n-ni, d-di) for ni, di in zip(ns,ds)]
//       fh = sum(fs)/l
//       vh = sum((f-fh)**2 for f in fs)/l
//       error2 = (l-1)*vh
//       return value, np.sqrt(error2)
//    else:
//       return value, 0

pub fn jacknife<F: Fn(f64, f64) -> f64>(f: F, aa: Vec<f64>, bb: Vec<f64>) -> (f64, f64) {
    let a: f64 = aa.iter().sum();
    let b: f64 = bb.iter().sum();
    let val = f(a, b);
    let l = aa.len();
    if l > 1 {
        let fs: Vec<f64> = aa
            .iter()
            .zip(bb.iter())
            .map(|(ia, ib)| f(a - ia, b - ib))
            .collect();
        let fh = fs.iter().cloned().sum::<f64>() / l as f64;
        let vh = fs.iter().map(|f| (f - fh).powi(2)).sum::<f64>() / l as f64;
        let err = (l - 1) as f64 * vh;
        (val, err)
    } else {
        (val, 0.0)
    }
}

pub fn vn(n: f64, events: &Vec<&Event>) -> (f64, f64) {
    let (nums, dens) = events
        .iter()
        .map(|e| {
            e.freezeouts
                .iter()
                .map(|f| {
                    let particles: Vec<_> = f
                        .particles
                        .iter()
                        .filter(|p| p.eta.abs() < 0.8 && 0.2 < p.pt && p.pt < 5.0)
                        .cloned()
                        .collect();

                    let m = particles.len() as f64;
                    let qn2 = particles
                        .iter()
                        .map(|p| Complex::from_polar(1.0, n * p.phi))
                        .sum::<Complex<f64>>()
                        .norm_sqr();

                    (qn2 - m, m * (m - 1.0))
                })
                .fold((0.0, 0.0), |(num, den), (m2, w2)| (num + m2, den + w2))
        })
        .fold((vec![], vec![]), |(mut nums, mut dens), (num, den)| {
            nums.push(num);
            dens.push(den);
            (nums, dens)
        });

    jacknife(|num, den| (num / den).sqrt(), nums, dens)
}

// Simetrized dn/deta os function of eta
pub fn dn_deta_eta(events: &Vec<&Event>) -> Vec<(f64, f64, f64)> {
    let deta = 0.25;
    let etas = (0..(5.0 / deta) as usize).map(|i| i as f64 * deta);

    etas.map(|eta| {
        let dndetas: Vec<f64> = events
            .iter()
            .map(|e| {
                e.freezeouts
                    .iter()
                    .map(|f| {
                        f.particles
                            .iter()
                            .filter(|p| eta <= p.eta.abs() && p.eta.abs() <= eta + deta)
                            .count() as f64
                    })
                    .fold(0.0, |acc, a| acc + a)
                    / e.freezeouts.len() as f64
                    / (2.0 * deta)
            })
            .collect();
        let dndeta: f64 = dndetas.iter().fold(0.0, |acc, a| acc + a) / dndetas.len() as f64;
        let dndeta2 = dndetas.iter().fold(0.0, |acc, a| acc + a * a) / dndetas.len() as f64;
        let errdndeta = ((dndeta2 - dndeta * dndeta) / events.len() as f64).sqrt();

        (eta + deta / 2.0, dndeta, errdndeta)
    })
    .collect()
}

pub fn dn_deta(events: &Vec<&Event>) -> (f64, f64) {
    let deta = 0.5;
    let dndetas: Vec<f64> = events
        .iter()
        .map(|e| {
            e.freezeouts
                .iter()
                .map(|f| f.particles.iter().filter(|p| p.eta.abs() <= deta).count() as f64)
                .fold(0.0, |acc, a| acc + a)
                / e.freezeouts.len() as f64
                / (2.0 * deta)
        })
        .collect();
    let dndeta: f64 = dndetas.iter().fold(0.0, |acc, a| acc + a) / dndetas.len() as f64;
    let dndeta2 = dndetas.iter().fold(0.0, |acc, a| acc + a * a) / dndetas.len() as f64;
    let errdndeta = ((dndeta2 - dndeta * dndeta) / events.len() as f64).sqrt();

    (dndeta, errdndeta)
}

fn main() {
    const DEBUG: bool = false;
    let mut events = vec![];
    for path in std::fs::read_dir("results").unwrap() {
        let p = path.unwrap().path();
        if DEBUG {
            println!("dir: {}", p.display());
        }

        let info: HashMap<String, String> = std::fs::read_to_string(p.join("info.txt"))
            .expect("Could not read info.txt")
            .trim()
            .split("\n")
            .map(|l| {
                let v: Vec<_> = l.split(": ").collect();
                (v[0].to_string(), v[1].to_string())
            })
            .collect();
        if DEBUG {
            println!("info: {:?}", info);
        }

        let freezeouts: Vec<FreezeOut> = std::fs::read_to_string(p.join("particles_out.dat"))
            .expect("Cloud not read particles_out.dat")
            .trim()
            .split("\n")
            .map(|l| l)
            .fold(vec![], |mut acc, l| {
                if l.starts_with("#") {
                    acc.push(FreezeOut { particles: vec![] });
                } else {
                    let last = acc.len() - 1;
                    // vals: ID charge pT ET mT phi y eta
                    //       0    1    2  3  4   5  6  7
                    let vals: Vec<f64> = l
                        .trim()
                        .split(" ")
                        .filter(|v| !v.is_empty())
                        .map(|v| v.parse().unwrap())
                        .collect();
                    if vals[1] != 0.0 {
                        // remove charless particles
                        acc[last].particles.push(Particle {
                            eta: vals[7],
                            pt: vals[2],
                            phi: vals[5],
                        });
                    }
                }
                acc
            });
        if DEBUG {
            println!("partirles: {:?}", freezeouts[0].particles[0]);
        }

        events.push(Event { freezeouts });

        if DEBUG {
            println!("");
        }
    }

    // Compute multiplicity of each events and attach the event index
    let mut mults: Vec<(usize, usize)> = events
        .iter()
        .enumerate()
        .map(|(i, e)| (e.multiplicity(), i))
        .collect();
    // Sort multiplicity to be able to bin. Larger multiplicity first
    mults.sort_by(|a, b| b.cmp(a));

    // Choose the number of bins
    let l = mults.len();
    let sizes = repeat(l / 40).take(4).chain(repeat(l / 10).take(9));
    let percent = |i| i as f64 / l as f64 * 100.0;

    let mut msg = vec![String::new(); 3];
    let mut j = 0;
    // For each bin
    for size in sizes {
        // Multiplicity range
        let bin = (percent(j), percent(j + size));

        // Take all the events corresponding to the current multiplicity bin
        let events: Vec<&Event> = (0..size).map(|k| &events[mults[j + k].1]).collect();

        let (dndeta, errdndeta) = dn_deta(&events);
        msg[0] = format!(
            "{}dn/deta-mid|{:.1}-{:.1}|{:.0}:{:.1}\n",
            msg[0], bin.0, bin.1, dndeta, errdndeta
        );
        let dndetaetas = dn_deta_eta(&events);
        let dndetaetas_data = dndetaetas
            .iter()
            .map(|(eta, dndeta, errdndeta)| format!("{:.3}:{:.0}:{:.1}", eta, dndeta, errdndeta))
            .collect::<Vec<_>>()
            .join(" ");
        msg[1] = format!(
            "{}dn/deta|{:.1}-{:.1}|{}\n",
            msg[1], bin.0, bin.1, dndetaetas_data
        );
        let (v2, errv2) = vn(2.0, &events);
        msg[2] = format!(
            "{}v2|{:.1}-{:.1}|{:.3}:{:.3}\n",
            msg[2], bin.0, bin.1, v2, errv2
        );

        j += size;
    }

    let msg = msg.join("\n");
    println!("{}", msg);
    std::fs::write("observables.txt", msg).expect("Could not write to observables.txt");
}
