use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
pub enum Integration {
    FixPoint,
    FixPointAnderson,
    Explicit,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Info {
    elapsed: f64,
    t0: f64,
    tend: f64,
    cost: u64,
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

    pub fn dim(&self) -> usize {
        let mut dim = 0;
        if self.nx > 1 {
            dim += 1;
        }
        if self.ny > 1 {
            dim += 1;
        }
        if self.nz > 1 {
            dim += 1;
        }
        dim
    }

    pub fn full_name(&self) -> String {
        let visc = match self.viscosity.as_str() {
            "Shear" => format!("Shear({})", self.etaovers),
            "Bulk" => format!("Bulk({})", self.zetaovers),
            "Both" => format!("Both({},{})", self.etaovers, self.zetaovers),
            a => format!("{a}"),
        };
        format!(
            "absolute_{}ref_{}D_{}_{}_{}t0_{}tend_{}dx_{}nx",
            self.scheme,
            self.dim(),
            self.name,
            visc,
            self.t0,
            self.tend,
            self.dx,
            self.nx
        )
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Analysis {
    name: String,
    dxs: Vec<Dxs>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Dxs {
    dx: f64,
    cases: Vec<Cases>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Cases {
    case: u64,
    integrations: Vec<Integrations>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Integrations {
    integration: Integration,
    analysis: Vec<Data>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Data {
    dt: f64,
    cost: u64,
    elapsed: f64,
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
            let info_str = std::fs::read_to_string(p.join(&stmax).join("info.txt")).unwrap();
            let info = serde_yaml::from_str::<Info>(&info_str)
                .unwrap()
                .post_process();
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

    let res = Analysis {
        name: values[0].0.full_name(),
        dxs: analysis(
            &values,
            |i| i.dx,
            |dx, _, values| Dxs {
                dx,
                cases: analysis(
                    values,
                    |i| i.case,
                    |case, _, values| Cases {
                        case,
                        integrations: analysis(
                            values,
                            |i| i.integration,
                            |integration, reference, values| Integrations {
                                integration,
                                analysis: {
                                    let data_ref = &reference[0].1;
                                    let e_tot_ref = data_ref.iter().sum::<f64>();
                                    values
                                        .iter()
                                        .map(|(info, d)| {
                                            let (etot, tot, max) =
                                                d.iter().zip(data_ref.iter()).fold(
                                                    (0.0f64, 0.0f64, 0.0f64),
                                                    |(etot, tot, max), (e, r)| {
                                                        let diff = (e - r).abs();
                                                        (etot + e, tot + diff, max.max(diff))
                                                    },
                                                );
                                            Data {
                                                dt: info.maxdt,
                                                elapsed: info.elapsed,
                                                cost: info.cost,
                                                relative_tot_energy: (e_tot_ref - etot).abs()
                                                    / e_tot_ref,
                                                mean_error: tot / d.len() as f64,
                                                max_error: max,
                                            }
                                        })
                                        .collect::<Vec<_>>()
                                },
                            },
                        ),
                    },
                ),
            },
        ),
    };
    let serialized = serde_yaml::to_string(&res).unwrap();
    std::fs::write("benchmark.txt", serialized).unwrap();
}

pub fn analysis<
    FO: PartialEq,
    Find: Fn(&Info) -> FO,
    D,
    Do: Fn(FO, &[(Info, Vec<f64>)], &[(Info, Vec<f64>)]) -> D,
>(
    mut values: &[(Info, Vec<f64>)],
    f: Find,
    d: Do,
) -> Vec<D> {
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
        res.push(d(start, first.unwrap(), &values[..pos]));
        values = &values[pos..];
    }
    res
}
