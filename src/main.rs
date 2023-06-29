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

    #[arg(short, long, default_value_t = Solver::Both)]
    solver_type: Solver,

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
        #[arg(long, default_value_t = 0.48)]
        t0: f64,
        #[arg(long, default_value_t = 10.0)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim3 {
        #[arg(long, default_value_t = 0.48)]
        t0: f64,
        #[arg(long, default_value_t = 3.0)]
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
        #[arg(long, default_value_t = 0.08)]
        etaovers_min: f64,
        #[arg(long, default_value_t = 1.23)] // GeV^-1
        etaovers_slope: f64,
        #[arg(long, default_value_t = -0.09)]
        etaovers_crv: f64,
        #[arg(long, default_value_t = 0.026)]
        zetaovers_max: f64,
        #[arg(long, default_value_t = 0.035)] // GeV
        zetaovers_width: f64,
        #[arg(long, default_value_t = 0.174)] // GeV
        zetaovers_peak: f64,
        #[arg(long, default_value_t = 0.050)] // GeV
        tempcut: f64,
        #[arg(long, default_value_t = 0.148)] // GeV
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
    let l = args.physical_length;
    let checkdt = |dtmin: f64, dtmax: f64| {
        if dtmin >= dtmax {
            panic!("The smallest time interval 'dtmin={}' must be smaller than the largest 'dtmax={}'.", dtmin, dtmax);
        }
    };
    let checkt = |t0: f64, tend: f64| {
        if t0 >= tend {
            panic!("The initial time 't0' must be smaller than the end time 'tend'.");
        }
    };
    let dx = l / XY as f64;
    let solver = args.solver_type;
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
                        let dtmax = dtmax.unwrap_or(dx * 0.1);
                        checkdt(dtmin, dtmax);
                        ideal1d::run_1d::<XY>(solver, t0, tend, l, dtmin, dtmax);
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
                        let dtmax = dtmax.unwrap_or(dx * 0.1);
                        checkdt(dtmin, dtmax);
                        ideal2d::run_2d::<XY>(
                            solver,
                            t0,
                            tend,
                            l,
                            dtmin,
                            dtmax,
                            enable_gubser,
                            nb_trento,
                        );
                    }
                    ToSimulate::Trento { dt, nb_trento } => {
                        let dt = dt.unwrap_or(dx * 0.1);
                        ideal2d::run_trento_2d::<XY>(solver, t0, tend, l, dt, nb_trento);
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
                            let dtmax = dtmax.unwrap_or(dx * 0.1);
                            checkdt(dtmin, dtmax);
                            viscous2d::run_2d::<XY>(
                                solver, t0, tend, l, dtmin, dtmax, etaovers, zetaovers, tempcut,
                                freezeout, nb_trento,
                            );
                        }
                        ToSimulate::Trento { dt, nb_trento } => {
                            let dt = dt.unwrap_or(dx * 0.1);
                            viscous2d::run_trento_2d::<XY>(
                                solver, t0, tend, l, dt, etaovers, zetaovers, tempcut, freezeout,
                                nb_trento,
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
                        let dtmax = dtmax.unwrap_or(dx * 0.1);
                        checkdt(dtmin, dtmax);
                        ideal3d::run_3d::<XY, Z>(solver, t0, tend, l, dtmin, dtmax, nb_trento);
                    }
                    ToSimulate::Trento { dt, nb_trento } => {
                        let dt = dt.unwrap_or(dx * 0.1);
                        ideal3d::run_trento_3d::<XY, Z>(solver, t0, tend, l, dt, nb_trento);
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
                            let dtmax = dtmax.unwrap_or(dx * 0.1);
                            checkdt(dtmin, dtmax);
                            viscous3d::run_3d::<XY, Z>(
                                solver, t0, tend, l, dtmin, dtmax, etaovers, zetaovers, tempcut,
                                freezeout, nb_trento,
                            );
                        }
                        ToSimulate::Trento { dt, nb_trento } => {
                            let dt = dt.unwrap_or(dx * 0.1);
                            viscous3d::run_trento_3d::<XY, Z>(
                                solver, t0, tend, l, dt, etaovers, zetaovers, tempcut, freezeout,
                                nb_trento,
                            );
                        }
                    }
                }
            }
        }
    }
}
