use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    let out_dir = format!("{}/bin", env::var("CARGO_HOME").unwrap());
    for name in ["implplt", "implPreHydro", "implPostHydro"] {
        let dest_path = PathBuf::from(&out_dir).join(name);
        println!("cargo:rerun-if-changed={name}.py");
        fs::copy(&format!("utils/{name}.py"), &dest_path)
            .expect(&format!("Failed to copy {name}.py."));
    }
}
