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

### Usefull scripts

Many useful scripts are located in the `utils` directory.  

_Description_:  
TODO