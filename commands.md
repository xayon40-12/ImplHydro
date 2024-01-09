TODO:
- use a general notation at the begining like: implhydro [general opt] {dim1|2|3} [dim opt] {Ideal|Viscous} [I|V opt] {Benchmark|Trento} [B|T opt]
- add the temperature dependent formula for $\eta$ and $\zeta$, and cite


# ImplHydro commands

This program use the convetion where a program can first have named arguments that start with a hyphen "-", and then one or many commands that are plain names without hyphen. Each command can recursively have named arguments and commands.  

### Notation:  
In this document, the program and commands will be surrounded by square brackets for readability. If a command can be followed by other sub-commands, the family of sub-commands (i.e. name of the section describing them) will specified in the title of the section of the command by separating the command and the sub-command family with a column ":". For instance, a section with title `[A] : B` will describe the command `[A]` which can/must use sub-commands described in the section `B`. Similarly, if some arguments uses special types of values that are not numeric, the unit will correspond to the name of a sub-section that will describe the avalable values.

When actually using the program, the square brackets should not be used. See the [Example](#Example) section.

---
## ImplHydro

### `[implhydro]` : [`Dim`](#Dim)
The main command.

| argument | description | default value | unit |
|----------|------|---------------|----|
| -c, --cells | Number of cells in the transverse X-Y plane. The number of cells in the $\eta$ direction is half this value. The 1D case is considered to be in the X direction, so no division by 2. | 100 | cells |
| -p, --physical_length | Total length of the system in the transverse X-Y plane. The length in the $\eta$ direction is half this value (although it is in unit of $\eta$, not in fm). The 1D case is considered to be in the X direction, so no division by 2. | 40 | fm |
| -e, --energy | Collision energy | 5020 | GeV |
| -s, --solver_type | Type of solver to be used. It can be explicit with `Explicit`, implicit with `FixPoint` or both with `Both` (this will generate one independent result for each case). | `Both` | [`Solver`](#Solver) |
| -r, --raw_data_every | Save the raw data every given interval in fm. If not used or used with the value `None`, no raw data are saved. | `None` | fm |

### `Solver`
| value | description |
|-------|-------------|
| `Explicit` | Explicit solver. |
| `FixPoint` | Implicit solver using a fixed-point iterator. |
| `Both` | Use both solver independently and solve the results in two different folders. |

---

## `Dim`
There are 3 dimensions available: 1D, 2D and 3D. For each dimension the initial and end time can be specified. If the end time in 0 fm, then the simulation will only end when all the cells are bellow the freeze-out temperature. Such freeze-out temperature can only be specified in the [`[Viscous]`](#Viscous) case of [`Hydro`](#Hydro), therefor the value 0 fm for the end time should only be used when simulating viscous hydrodynamics.

### `[Dim1]` : [`Hydro`](#Hydro)
1D case.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --t0 | Initial time | 0 | fm |
| --tend | End time | 15 | fm |

### `[Dim2]` : [`Hydro`](#Hydro)
2D case.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --t0 | Initial time | 0.37 | fm |
| --tend | End time | 0 | fm |

### `[Dim3]` : [`Hydro`](#Hydro)
3D case.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --t0 | Initial time | 0.37 | fm |
| --tend | End time | 0 | fm |

---

## `Hydro`
Choose whether to use ideal of viscous hydrodynamics.

### `[Ideal]` : [`ToSimulate`](#ToSimulate)
Ideal hydrodynamics.

| argument | description | default value | unit |
|----------|------|---------------|----|
| -e, --enable_gubser | In addition to the TrENTo initial conditions (which can be disabled in the [`ToSimulate`](#ToSimulate) sub-commands), the Gubser flow solution is used as initial conditions for a simulation that will be performed and saved independently of the TrENTo one. | false | true/false |

### `[Viscous]` : [`ToSimulate`](#ToSimulate)
Viscous hydrodynamics.  
The viscous transport coefficients $\eta$ and $\zeta$ are temperature dependent.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --etaovers_min |   | 0.11 |   |
| --etaovers_slope |   | 1.6 | 1/GeV |
| --etaovers_crv |   | -0.29 |   |
| --zetaovers_max |   | 0.032 |   |
| --zetaovers_width |   | 0.024 | GeV |
| --zetaovers_peak |   | 0.175 | GeV |
| --tempcut | Temperature used to stabilize the numerical hydrodynamics. Below this temperature, the viscosity is disabled and relax to 0 instead of its Navier-Stocks value. Without such regulation, the viscosity would diverge near the vacuum. | 0.020 | GeV |

---

## `ToSimulate`
Choose whether to perform a time accuracy benchmark by varying $\Delta t$ or just to run for a single $Delta t$.

### `[Benchmark]`
Run many simulations by dividing $Delta t$ by 2 many times, starting at `dtmax` until it reaches `dtmin`.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --dtmin | Minimum value of $\Delta t$. As $\Delta t$ is divided many times by 2 starting from `dtmax`, the value `dtmin` might not be reached exactly. | 1e-3 | fm |
| --dtmax | Maximum value of $\Delta t$ and initial value | 0.2*physical_length/cells | fm |
| -n, --nb_trento | Number of TrENTo initial conditions to use. | 10 |   |

### `[Trento]`
Run the simulation for a fixed $\Delta t$.

| argument | description | default value | unit |
|----------|------|---------------|----|
| --dt | $\Delta t$ | 0.2*physical_length/cells | fm |
| -n, --nb_trento | Number of TrENTo initial conditions to use. | 10 |   |

## Example
To run a (1+1)-dimensional benchmark fo the Riemann problem by saving the data every fm of time:  
```bash
implhydro -c 200 -r 1 dim1 ideal benchmark --dtmin 5e-4 --dtmax 1.28
```
