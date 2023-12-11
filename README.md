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
| `implPreHydro.py` | Call [`trento`](https://github.com/Duke-QCD/trento) or [`trento3d-2.0`](https://github.com/Duke-QCD/trento3d-2.0) to generate the initial condition. The flag `-3d` can be used to enable 3d `trento`. The `-f` flag can be used to use the [`freestream`](https://github.com/Duke-QCD/freestream) pre-hydrodynamics. |
| `implPostHydro.py` | Call [`frzout`](https://github.com/Duke-QCD/frzout) on the `surface.dat` hyper-surface output (must be called directly next to such file in the `results`), and then call the [`UrQMD`](https://github.com/jbernhard/urqmd-afterburner) afterburner. |
| `plt_setting.py` | Settings for matplotlib from [Saizensen](https://github.com/MasakiyoK/Saizensen) |

## File formats

### Input format

The `implhydro` code expect the input data to be formated in the same way that `trento` output its readable format (hdf5 is not supported). The data files should be stored in a folder that starts by the letter 's' followed by the number of cells in the X-Y directions, where the X and Y directions must have the same number of cells. For instance, for a 100 by 100 cells lattice, the data files should be stored in a folder called `s100`. The data files should be named by a number with left zero padding followed by the `.dat` file extension. For instance, in the case where 100 initial conditions are stored, the first file would be `00.dat` and the last one `99.dat`.  

#### 2+1D matrix format

For the 2+1D case, the columns correpsond to the X direction and the lines to the Y direction, where the X-Y plane correspond to the collision plane or transverse plane in a heavy-ion collision.

#### 3+1D matrix format

For the 3+1D case, the columns correspond to the X direction and the lines to the Y and Z/$\eta$ direction, where the X-Y plane correspond to the collision plane or tranverse plone in a heavy-ion collision. The Y direction comes first and is repeated for each cells in the Z/$\eta$ direction:  

- line: Z, Y
-   01: 0, 0
-   02: 0, 1
- ...
-   10: 0, 9
-   11: 1, 0
-   12: 1, 2
- ...