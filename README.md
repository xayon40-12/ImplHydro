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

Then, build and install the executable `implhydro`:  
```sh
cargo install --path .
```

### Usefull scripts

Many useful scripts are located in the `utils` directory. The most useful ones are installed as well without the `.py` extension during the install process of the `implhydro` executable for easier access.
