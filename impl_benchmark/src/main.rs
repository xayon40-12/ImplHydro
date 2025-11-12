use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
pub enum Integration {
    FixPoint,
    FixPointAnderson,
    Explicit,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Info {
    t0: f64,
    tend: f64,
    nx: u64,
    ny: u64,
    nz: u64,
    dx: f64,
    dy: f64,
    dz: f64,
    maxdt: f64,
    integration: Integration,
    scheme: String,
    stages: u64,
    name: String,
    case: u64,
    viscosity: String,
    etaovers: String,
    zetaovers: String,
    energycut: f64,
    freezeoutenergy: f64,
    variables: String,
    nb_variables: Option<usize>,
    energy_index: Option<usize>,
}

impl Info {
    pub fn post_process(mut self) -> Self {
        self.energy_index = self.variables.split(" ").position(|x| x == "e");
        self.nb_variables = Some(self.variables.chars().filter(|&c| c == ' ').count() + 1);
        self
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Analysis {
    dt: f64,
    relative_tot_energy: f64,
    mean_error: f64,
    max_error: f64,
}

fn main() {
    let mut values = std::fs::read_dir("results")
        .unwrap()
        .into_iter()
        .map(|path| {
            let p = path.unwrap().path();
            let info_str = std::fs::read_to_string(p.join("info.txt")).unwrap();
            let info = serde_yaml::from_str::<Info>(&info_str)
                .unwrap()
                .post_process();
            let stmax = std::fs::read_dir(&p)
                .unwrap()
                .filter_map(|d| {
                    d.ok()
                        .and_then(|d| d.file_name().to_str().map(|s| s.to_string()))
                        .and_then(|d| d.parse::<f32>().ok().map(|x| (x, d)))
                })
                .max_by(|(x, _), (y, _)| x.total_cmp(y))
                .unwrap()
                .1;
            let p = p.join(&stmax).join("data.dat");

            let eid = info.energy_index.unwrap();
            let data = std::fs::read(&p)
                .unwrap()
                .chunks(8)
                .map(|cs| unsafe { std::mem::transmute::<&[u8], &[f64]>(cs) }[0])
                .collect::<Vec<f64>>()
                .chunks(info.nb_variables.unwrap())
                .map(|cs| cs[eid])
                .collect::<Vec<f64>>();
            (info, data)
        })
        .collect::<Vec<_>>();

    values.sort_by(|(a, _), (b, _)| {
        a.dx.total_cmp(&b.dx)
            .reverse()
            .then(a.case.cmp(&b.case))
            .then(a.integration.cmp(&b.integration))
            .then(a.maxdt.total_cmp(&b.maxdt))
    });

    let res = analysis(
        &values,
        |i| i.dx,
        |_, values| {
            analysis(
                values,
                |i| i.case,
                |_, values| {
                    analysis(
                        values,
                        |i| i.integration,
                        |reference, values| {
                            let data_ref = &reference[0].1;
                            let e_tot_ref = data_ref.iter().sum::<f64>();
                            values
                                .iter()
                                .map(|(info, d)| {
                                    let (etot, tot, max) = d.iter().zip(data_ref.iter()).fold(
                                        (0.0f64, 0.0f64, 0.0f64),
                                        |(etot, tot, max), (e, r)| {
                                            let diff = (e - r).abs();
                                            (etot + e, tot + diff, max.max(diff))
                                        },
                                    );
                                    Analysis {
                                        dt: info.maxdt,
                                        relative_tot_energy: (e_tot_ref - etot).abs() / e_tot_ref,
                                        mean_error: tot / d.len() as f64,
                                        max_error: max,
                                    }
                                })
                                .collect::<Vec<_>>()
                        },
                    )
                },
            )
        },
    );
    println!("{res:?}");
}

pub fn analysis<
    FO: PartialEq,
    Find: Fn(&Info) -> FO,
    D,
    Do: Fn(&[(Info, Vec<f64>)], &[(Info, Vec<f64>)]) -> D,
>(
    mut values: &[(Info, Vec<f64>)],
    f: Find,
    d: Do,
) -> Vec<(FO, D)> {
    let mut res = vec![];
    let mut first = None;
    while values.len() > 0 {
        let start = f(&values[0].0);
        let pos = values
            .iter()
            .position(|v| f(&v.0) != start)
            .unwrap_or(values.len());
        if first.is_none() {
            first = Some(&values[..pos]);
        }
        res.push((start, d(first.unwrap(), &values[..pos])));
        values = &values[pos..];
    }
    res
}
