use clap::{command, Parser, Subcommand};
use implhydro::{
    run::{ideal1d, ideal2d, ideal3d, viscous2d, viscous3d},
    solver::Solver,
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(short, long, default_value_t = 100)]
    cells: usize,

    #[arg(short, long, default_value_t = 40.0)]
    physical_length: f64,

    #[arg(short, long, default_value_t = 5020.0)] // GeV
    collision_energy: f64,

    #[arg(short, long, default_value_t = Solver::Both)]
    solver_type: Solver,

    #[arg(short, long)]
    raw_data_every: Option<f64>,

    #[command(subcommand)]
    command: Dim,
}

#[derive(Subcommand, Debug)]
enum Dim {
    Dim1 {
        #[arg(long, default_value_t = 0.0)]
        t0: f64,
        #[arg(long, default_value_t = 15.0)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim2 {
        #[arg(long, default_value_t = 0.37)]
        t0: f64,
        #[arg(long, default_value_t = 0.0)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim3 {
        #[arg(long, default_value_t = 1.3)]
        t0: f64,
        #[arg(long, default_value_t = 0.0)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
}

#[derive(Subcommand, Debug)]
enum Hydro {
    Ideal {
        #[arg(short, long, default_value_t = false)]
        enable_gubser: bool,
        #[command(subcommand)]
        command: ToSimulate,
    },
    Viscous {
        #[arg(long, default_value_t = 0.11)]
        etaovers_min: f64,
        #[arg(long, default_value_t = 1.6)] // GeV^-1
        etaovers_slope: f64,
        #[arg(long, default_value_t = -0.29)]
        etaovers_crv: f64,
        #[arg(long, default_value_t = 0.032)]
        zetaovers_max: f64,
        #[arg(long, default_value_t = 0.024)] // GeV
        zetaovers_width: f64,
        #[arg(long, default_value_t = 0.175)] // GeV
        zetaovers_peak: f64,
        #[arg(long, default_value_t = 0.020)] // GeV
        tempcut: f64,
        #[arg(long, default_value_t = 0.151)] // GeV
        freezeout: f64,
        #[command(subcommand)]
        command: ToSimulate,
    },
}

#[derive(Subcommand, Debug)]
enum ToSimulate {
    Benchmark {
        #[arg(long, default_value_t = 1e-3)]
        dtmin: f64,
        #[arg(long)]
        dtmax: Option<f64>,
        #[arg(short, long, default_value_t = 10)]
        nb_trento: usize,
    },
    Trento {
        #[arg(long)]
        dt: Option<f64>,
        #[arg(short, long, default_value_t = 100)]
        nb_trento: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.cells {
        100 => run::<100, 51>(cli),
        200 => run::<200, 101>(cli),
        _ => panic!("The number of cells must be a value from the list {{100,200}}."),
    };
}
fn run<const XY: usize, const Z: usize>(args: Cli) {
    let sqrts = args.collision_energy;
    let xy_len = args.physical_length;
    let etas_len = 2.0 * (0.5 * sqrts / 0.2).acosh();
    let checkdt = |dtmin: f64, dtmax: f64| {
        if dtmin >= dtmax {
            panic!("The smallest time interval 'dtmin={}' must be smaller than the largest 'dtmax={}'.", dtmin, dtmax);
        }
    };
    let checkt = |_t0: f64, _tend: f64| {
        // if t0 > tend {
        //     panic!("The initial time 't0' must be smaller than the end time 'tend'.");
        // }
    };
    let dx = xy_len / XY as f64;
    let dtdx = 0.2f64;
    let solver = args.solver_type;
    let save_raw = args.raw_data_every;
    match args.command {
        Dim::Dim1 { t0, tend, command } => {
            checkt(t0, tend);
            match command {
                Hydro::Ideal {
                    enable_gubser: _,
                    command,
                } => match command {
                    ToSimulate::Benchmark {
                        dtmin,
                        dtmax,
                        nb_trento: _,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dtdx);
                        checkdt(dtmin, dtmax);
                        ideal1d::run_1d::<XY>(solver, t0, tend, xy_len, dtmin, dtmax, save_raw);
                    }
                    ToSimulate::Trento { .. } => panic!("No 1D Trento case available."),
                },
                Hydro::Viscous { .. } => match command {
                    _ => panic!("No 1D viscous case available."),
                },
            }
        }
        Dim::Dim2 { t0, tend, command } => {
            checkt(t0, tend);
            match command {
                Hydro::Ideal {
                    enable_gubser,
                    command,
                } => match command {
                    ToSimulate::Benchmark {
                        dtmin,
                        dtmax,
                        nb_trento,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dtdx);
                        checkdt(dtmin, dtmax);
                        ideal2d::run_2d::<XY>(
                            solver,
                            t0,
                            tend,
                            xy_len,
                            dtmin,
                            dtmax,
                            enable_gubser,
                            nb_trento,
                            save_raw,
                        );
                    }
                    ToSimulate::Trento { dt, nb_trento } => {
                        let dt = dt.unwrap_or(dx * dtdx);
                        ideal2d::run_trento_2d::<XY>(
                            solver, t0, tend, xy_len, dt, nb_trento, save_raw,
                        );
                    }
                },
                Hydro::Viscous {
                    etaovers_min,
                    etaovers_slope,
                    etaovers_crv,
                    zetaovers_max,
                    zetaovers_width,
                    zetaovers_peak,
                    tempcut,
                    freezeout,
                    command,
                } => {
                    let etaovers = (etaovers_min, etaovers_slope, etaovers_crv);
                    let zetaovers = (zetaovers_max, zetaovers_width, zetaovers_peak);
                    match command {
                        ToSimulate::Benchmark {
                            dtmin,
                            dtmax,
                            nb_trento,
                        } => {
                            let dtmax = dtmax.unwrap_or(dx * dtdx);
                            checkdt(dtmin, dtmax);
                            viscous2d::run_2d::<XY>(
                                solver, t0, tend, xy_len, dtmin, dtmax, etaovers, zetaovers,
                                tempcut, freezeout, nb_trento, save_raw,
                            );
                        }
                        ToSimulate::Trento { dt, nb_trento } => {
                            let dt = dt.unwrap_or(dx * dtdx);
                            viscous2d::run_trento_2d::<XY>(
                                solver, t0, tend, xy_len, dt, etaovers, zetaovers, tempcut,
                                freezeout, nb_trento, save_raw,
                            );
                        }
                    }
                }
            }
        }
        Dim::Dim3 { t0, tend, command } => {
            checkt(t0, tend);
            match command {
                Hydro::Ideal {
                    enable_gubser: _,
                    command,
                } => match command {
                    ToSimulate::Benchmark {
                        dtmin,
                        dtmax,
                        nb_trento,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dtdx);
                        checkdt(dtmin, dtmax);
                        ideal3d::run_3d::<XY, Z>(
                            solver, t0, tend, xy_len, etas_len, dtmin, dtmax, nb_trento, save_raw,
                        );
                    }
                    ToSimulate::Trento { dt, nb_trento } => {
                        let dt = dt.unwrap_or(dx * dtdx);
                        ideal3d::run_trento_3d::<XY, Z>(
                            solver, t0, tend, xy_len, etas_len, dt, nb_trento, save_raw,
                        );
                    }
                },
                Hydro::Viscous {
                    etaovers_min,
                    etaovers_slope,
                    etaovers_crv,
                    zetaovers_max,
                    zetaovers_width,
                    zetaovers_peak,
                    tempcut,
                    freezeout,
                    command,
                } => {
                    let etaovers = (etaovers_min, etaovers_slope, etaovers_crv);
                    let zetaovers = (zetaovers_max, zetaovers_width, zetaovers_peak);
                    match command {
                        ToSimulate::Benchmark {
                            dtmin,
                            dtmax,
                            nb_trento,
                        } => {
                            let dtmax = dtmax.unwrap_or(dx * dtdx);
                            checkdt(dtmin, dtmax);
                            viscous3d::run_3d::<XY, Z>(
                                solver, t0, tend, xy_len, etas_len, dtmin, dtmax, etaovers,
                                zetaovers, tempcut, freezeout, nb_trento, save_raw,
                            );
                        }
                        ToSimulate::Trento { dt, nb_trento } => {
                            let dt = dt.unwrap_or(dx * dtdx);
                            viscous3d::run_trento_3d::<XY, Z>(
                                solver, t0, tend, xy_len, etas_len, dt, etaovers, zetaovers,
                                tempcut, freezeout, nb_trento, save_raw,
                            );
                        }
                    }
                }
            }
        }
    }
}
