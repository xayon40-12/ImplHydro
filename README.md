# ImplHydro

ImplHydro is a 1+1D, 2+1D and 3+1D implicit in time relativistic hydrodynamics solver. It also includes an explicit solver for completeness.  
Its purpose is to compare the efficiency of explicit and implicit Runge-Kutta 
method used to integrate the conservation equation of the energy momentum tensor $T^{\mu\nu}$.

## Installation


### Rust compiler
In order to execute this project, you will need a working `rust` compiler along with its project manager called `cargo`.  
The simplest way to have both is to use the official installer called `rustup` available at [rustup.rs](https://rustup.rs).  


### implhydro
Once `rust` is installed, you can clone this `ImplHydro` project and go inside its directory:  
```bash
git clone https://github.com/xayon40-12/ImplHydro.git
cd ImplHydro
```

If you installed `rust` with `rustup`, then you should have a `.cargo` directory in your home folder. Furthermore, you should have been asked to add `$HOME/.cargo/bin` to you `PATH` environement variable. If so, you can install the `implhydro` executable and the other useful scripts located in the `utils` folder by executing the `install` script:  
```bash
./install
```

If you do not want to install, you can simply build the project with
```bash
cargo build --release
```
and then acces the executable in `target/release/implhydro`

## Usage

### Run simulations

Once the `implhydro` executable is compiled and available, you can simply call it with the `--help` flag
```bash
implhydro --help
```
and it will show a help message that explain the available options with their default values and the different commands that are to be used.
Every command has it own options and commands. The same `--help` flag can be used to obtain explanations for every sub-commands.  

To reproduce the 1-dimensional case results and then generate the figures, you can run:
```bash
implhydro -c 200 -r 1 dim1 ideal benchmark --dtmin 5e-4 --dtmax 1.28
implhydro -c 100 -r 1 dim1 ideal benchmark --dtmin 5e-4 --dtmax 1.28
implplt.py -r
```

TODO:  
- detail the subcommands as the default values do not detail what is available (and some defaults are not cli defaults).
- give an exemple of Trento3D

### Usefull scripts

Many useful scripts are located in the `utils` directory.  

| Name | Description |
|------|-------------|
| `implplt.py` | If executed in the same folder as the `results` folder, creates a `figures` folder and generates all the figures. The `-r` flag can be used to remove failed points, the `-s` flag shows scheme names, the `-a` flag save a time animation, the `-m` flag save figures for all the Trento results and not only the one labeled `0`. |
| `implPreHydro.py` | Call [`trento`](https://github.com/Duke-QCD/trento) or [`trento3d-2.0`](https://github.com/Duke-QCD/trento3d-2.0) to generate the initial condition. The flag `-3d` can be used to enable 3d `trento`. |
| `implPostHydro.py` | Must be used next to the `results` folder. It will execute [`frzout`](https://github.com/Duke-QCD/frzout) on all the `surface.dat` hyper-surface outputs to create a `particles_in.dat` file, then it will execute the [`UrQMD`](https://github.com/jbernhard/urqmd-afterburner) afterburner to create the `particles_out.dat` file which contains the final particles.|
| `plt_setting.py` | Settings for matplotlib from [Saizensen](https://github.com/MasakiyoK/Saizensen) |

In addition, there is another rust program that is install with the `./install` command called `analysis`. After running the `implPostHydro.py` script, the `analysis` progam can be executed next to the `results` folder to generate a file named `observables.txt` which will contains the observables $\mathrm{dN}_\mathrm{ch}/\mathrm{d}\eta$ and $v_2\{2\}$ as function of centrality.  

## File formats

### Input format

The `implhydro` code expect the input data to be formated in the same way that `trento` outputs its readable format (hdf5 is not supported). The data files should be stored in a folder that starts by the letter 's' followed by the number of cells in the X-Y directions, where the X and Y directions must have the same number of cells. For instance, for a 100 by 100 cells lattice, the data files should be stored in a folder called `s100`. The data files should be named by a number with left zero padding followed by the `.dat` file extension. For instance, in the case where 100 initial conditions are stored, the first file would be `00.dat` and the last one `99.dat`.  

#### 2+1D matrix format

For the 2+1D case, the columns correpsond to the X direction and the lines to the Y direction, where the X-Y plane correspond to the collision plane or transverse plane in a heavy-ion collision.

#### 3+1D matrix format

For the 3+1D case, the columns correspond to the X direction and the lines to the Y and Z/$\eta$ direction, where the X-Y plane correspond to the collision plane or tranverse plone in a heavy-ion collision. The Y direction comes first and is repeated for each cells in the Z/$\eta$ direction.

### Output format

The generated data are outputed in a folder named after the case simulated and the main parameters values. This folder is itself stored in the `result` folder.  
A file named `info.txt` is stored next to the outputed data. It contains all the informations involved such as initial time, lattice spacing, ...

#### Raw data

The raw data can be outputed by using the `-r <every>` flag where `<every>` is the physical time between each save. The raw data is a binary file of double precison floating point numbers.   
All the data of one cell are stored together and every cell come in order where the X direction is first, Y second and Z/$\eta$ last. The `variables:` field of the `info.txt` file contains the name of the variables stored in every cells. There will be one double precision floating point number for each of these variables. As the first variables for each cells are the coordinates x, y, z, it is not necessary to take care about the order in which each cells are stored

#### Hypersurface data

In the case of 2+1D and 3+1D viscous with TrENTO initial conditions, the hypersurface data needed for the freezeout are stored in binary format in the `surface.dat` file. The data are stored in order as double precision floating point numbers: time $\tau$, space coordinate $x-y-\eta$, covariant normal vector $\sigma_\mu$ with length equal to the transverse volume, velocity $v^i$,  upper triangle part of the shear tensor$\pi^{\mu\nu}$, shear pressure.  
In 2+1D it corresponds to:  


| $\tau$ | $x$ | $y$ | $\sigma_\tau$ | $\sigma_x$ | $\sigma_y$ | $v^x$ | $v^y$ | $\pi^{\tau\tau}$ | $\pi^{\tau x}$ | $\pi^{\tau y}$ | $\pi^{x x}$ | $\pi^{xy}$ | $\pi^{yy}$ | $\pi^{\eta\eta}$ | $\Pi$ |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|

In 3+1D it correpsonds to:  

| $\tau$ | $x$ | $y$ | $\eta$ | $\sigma_\tau$ | $\sigma_x$ | $\sigma_y$ | $\sigma_\eta$ | $v^x$ | $v^y$ | $v^\eta$ | $\pi^{\tau\tau}$ | $\pi^{\tau x}$ | $\pi^{\tau y}$ | $\pi^{\tau\eta}$ | $\pi^{xx}$ | $\pi^{xy}$ | $\pi^{x\eta}$ | $\pi^{yy}$ | $\pi^{y\eta}$ | $\pi^{\eta\eta}$ | $\Pi$ |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
