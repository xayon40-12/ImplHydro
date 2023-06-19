# ImplHydro

ImplHydro is a 1+1D and 2+1D implicit in time relativistic hydrodynamics solver.  
Its purpose is to compare the efficiency of explicit and implicit Runge-Kutta 
method used to integrate the conservation equation of the energy momentum tensor $T^{\mu\nu}$.

## Installation


### Rust compiler
In order to execute this project, you will need a working `rust` compiler along with its project manager called `cargo`.  
The simplest way to have both is to use the official installer called `rustup` available at [rustup.rs](https://rustup.rs).  


### implhydro
Once `rust` is installed, you can clone this `ImplHydro` project and go inside its directory:  
```sh
git clone https://github.com/xayon40-12/ImplHydro.git
cd ImplHydro
```

If you installed `rust` with `rustup`, then you should have a `.cargo` directory in your home folder. Furthermore, you should have been asked to add `$HOME/.cargo/bin` to you `PATH` environement variable. If so, you can install the `implhydro` executable and the other useful scripts located in the `utils` folder by executing the `install` script:  
```sh
./install
```

If you do not want to install, you can simply build the project with
```sh
cargo build --release
```
and then acces the executable in `target/release/implhydro`

## Usage

### Run simulations

Once the `implhydro` executable is compiled and available, you can simply call it with the `--help` flag
```sh
implhydro --help
```
and it will show a help message that explain the available options with their default values and the different commands that are to be used.
Every command has it own options and commands. The same `--help` flag can be used to obtain explanations for every sub-commands.  

To reproduce the 1-dimensional case results and then generate the figures, you can run:
```sh
implhydro -c 200 dim1 ideal benchmark --dtmin 5e-4 --dtmax 1.28
implhydro -c 100 dim1 ideal benchmark --dtmin 5e-4 --dtmax 1.28
implplt.py -r
```

### Usefull scripts

Many useful scripts are located in the `utils` directory.  

| Name | Description |
|------|-------------|
| `implplt.py` | If executed in the same folder as the `results` folder, creates a `figures` folder and generates all the figures. The `-r` flag can be used to remove failed points, the `-s` flag show scheme names, the `-a` flag save a time animation, the `-m` flag save figures for all the Trento results and not only the one labeled `0`. |
| `implPreHydro.py` | Call [`TrENTo`](https://github.com/Duke-QCD/trento) to generate the initial condition. The `-f` flag can be used to use the [`freestream`](https://github.com/Duke-QCD/freestream) pre-hydrodynamics. |
| `implPostHydro.py` | Call [`frzout`](https://github.com/Duke-QCD/frzout) on the `surface.dat` hyper-surface output (must be called directly next to such file in the `results`), and then call the [`UrQMD`](https://github.com/jbernhard/urqmd-afterburner) afterburner. |
