use serde::{Deserialize, Serialize};

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
    integration: String,
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
    energy_index: Option<usize>,
}

impl Info {
    pub fn post_process(mut self) -> Self {
        self.energy_index = self.variables.split(" ").position(|x| x == "e");
        self
    }
}

fn main() {
    for path in std::fs::read_dir("results").unwrap() {
        let p = path.unwrap().path();
        let info_str = std::fs::read_to_string(p.join("info.txt")).unwrap();
        let info = serde_yaml::from_str::<Info>(&info_str)
            .unwrap()
            .post_process();
        println!("{info:?}");
        for times in std::fs::read_dir(p).unwrap().filter_map(|d| {
            d.ok()
                .and_then(|d| d.file_name().to_str().map(|s| s.to_string()))
                .and_then(|d| d.parse::<f32>().ok().map(|_| d))
        }) {
            println!("{times}");
        }
    }
}
