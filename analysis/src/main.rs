use std::{collections::HashMap, iter::repeat};

#[derive(Debug)]
pub struct Particle {
    pub eta: f64,
    pub pt: f64,
    pub phi: f64,
}

#[derive(Debug)]
pub struct FreezeOut {
    pub particles: Vec<Particle>,
}

#[derive(Debug)]
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

// def vn(n, events):
//    events = [[[v for v in f if abs(v["eta"])<0.8 and 0.2 < v["pT"] and v["pT"] < 5.0 and v["charge"] != 0] for f in e] for e in events]

//    def fun(n, d):
//       return np.sqrt(n/d)

//    mss = [[len(f) for f in e] for e in events]

//    qnss = [[sum(exp(1j*n*v["phi"]) for v in f) for f in e] for e in events]
//    qn2ss = [[abs(qn)**2 for qn in qns] for qns in qnss]

//    w2ss = [[m*(m-1) for m in ms] for ms in mss]

//    m2ss = [[(qn2-m)/w2 for qn2, m, w2 in zip(qn2s, ms, w2s)] for qn2s, ms, w2s in zip(qn2ss, mss, w2ss)]

//    n_cn2s = [sum(w2*m2 for w2, m2 in zip(w2s, m2s)) for w2s, m2s in zip(w2ss, m2ss)]
//    d_cn2s = [sum(w2s) for w2s in w2ss]

//    val, err = jackknife(fun, n_cn2s, d_cn2s)

//    return val, err

// def dch_deta_eta(events):
//    deta = 0.25
//    ran = np.arange(0, 5, deta)
//    etas = [eta+deta/2 for eta in ran]
//    dchdetass = [[len([v for v in f if eta <= abs(v["eta"]) and abs(v["eta"]) <= eta+deta and v["charge"] != 0])/(2*deta) for e in events for f in e] for eta in ran] # abs and /2 because it is symetrized
//    dchdetas = [sum(dchdetas)/len(dchdetas) for dchdetas in dchdetass]
//    dchdetas2 = [sum(x**2 for x in dchdetas)/len(dchdetas) for dchdetas in dchdetass]
//    errdchdetas = [np.sqrt((d2-d**2)/len(events)) for d, d2 in zip(dchdetas,dchdetas2)]
//    plt.errorbar(etas, dchdetas, yerr=errdchdetas)
//    plt.ylim([0,2500])
//    plt.savefig("dch_deta_eta.pdf")

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
    let sizes = repeat(l / 40).take(4).chain(repeat(l / 10).take(7));
    let percent = |i| i as f64 / l as f64 * 100.0;

    let mut msg = String::new();
    let mut j = 0;
    // For each bin
    for size in sizes {
        // Multiplicity range
        let bin = (percent(j), percent(j + size));

        // Take all the events corresponding to the current multiplicity bin
        let events: Vec<&Event> = (0..size).map(|k| &events[mults[j + k].1]).collect();

        let (dndeta, errdndeta) = dn_deta(&events);
        msg = format!(
            "{}dn/deta-mid|{:.1}-{:.1}|{:.0}:{:.1}\n",
            msg, bin.0, bin.1, dndeta, errdndeta
        );
        let dndetaetas = dn_deta_eta(&events);
        let dndetaetas_data = dndetaetas
            .iter()
            .map(|(eta, dndeta, errdndeta)| format!("{:.3}:{:.0}:{:.1}", eta, dndeta, errdndeta))
            .collect::<Vec<_>>()
            .join(" ");
        msg = format!(
            "{}dn/deta|{:.1}-{:.1}|{}\n",
            msg, bin.0, bin.1, dndetaetas_data
        );

        j += size;
    }
    println!("{}", msg);
    std::fs::write("observables.txt", msg).expect("Could not write to observables.txt");
}
