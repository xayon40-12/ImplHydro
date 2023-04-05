use clap::{command, Parser, Subcommand};
use implhydro::run::{ideal, shear};
use std::thread;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(short, long, default_value_t = 100)]
    cells: usize,

    #[arg(short, long, default_value_t = 20.0)]
    physical_length: f64,

    #[command(subcommand)]
    command: Dim,
}

#[derive(Subcommand, Debug)]
enum Dim {
    Dim1 {
        #[arg(long, default_value_t = 0.0)]
        t0: f64,
        #[arg(long, default_value_t = 7.5)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim2 {
        #[arg(long, default_value_t = 0.48)]
        t0: f64,
        #[arg(long, default_value_t = 15.0)]
        tend: f64,
        #[command(subcommand)]
        command: Hydro,
    },
}

#[derive(Subcommand, Debug)]
enum Hydro {
    Ideal {
        #[command(subcommand)]
        command: ToSimulate,
    },
    Shear {
        #[arg(short, long, default_value_t = 0.08)]
        etaovers_min: f64,
        #[arg(short, long, default_value_t = 1.23e-3)] // MeV^-1
        etaovers_slope: f64,
        #[arg(short, long, default_value_t = -0.09)]
        etaovers_crv: f64,
        #[arg(short, long, default_value_t = 0.026)]
        zetaovers_max: f64,
        #[arg(short, long, default_value_t = 35.0)] // MeV
        zetaovers_width: f64,
        #[arg(short, long, default_value_t = 174.0)] // MeV
        zetaovers_peak: f64,
        #[arg(long, default_value_t = 20.0)]
        tempcut: f64,
        #[arg(long, default_value_t = 148.0)]
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

fn big_stack() {
    let cli = Cli::parse();

    match cli.cells {
        100 => run::<100>(cli),
        200 => run::<200>(cli),
        _ => panic!("The number of cells must be a value from the list {{100,200}}."),
    };
}
fn run<const CELLS: usize>(args: Cli) {
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
    let dx = l / CELLS as f64;
    match args.command {
        Dim::Dim1 { t0, tend, command } => {
            checkt(t0, tend);
            match command {
                Hydro::Ideal { command } => match command {
                    ToSimulate::Benchmark {
                        dtmin,
                        dtmax,
                        nb_trento: _,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * 0.1);
                        checkdt(dtmin, dtmax);
                        ideal::run_1d::<CELLS>(t0, tend, l, dtmin, dtmax);
                    }
                    ToSimulate::Trento { .. } => panic!("No 1D Trento case available."),
                },
                Hydro::Shear { .. } => match command {
                    _ => panic!("No 1D shear case available."),
                },
            }
        }
        Dim::Dim2 { t0, tend, command } => {
            checkt(t0, tend);
            match command {
                Hydro::Ideal { command } => match command {
                    ToSimulate::Benchmark {
                        dtmin,
                        dtmax,
                        nb_trento,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * 0.1);
                        checkdt(dtmin, dtmax);
                        ideal::run_2d::<CELLS>(t0, tend, l, dtmin, dtmax, nb_trento);
                    }
                    ToSimulate::Trento { dt, nb_trento } => {
                        let dt = dt.unwrap_or(dx * 0.1);
                        ideal::run_trento_2d::<CELLS>(t0, tend, l, dt, nb_trento);
                    }
                },
                Hydro::Shear {
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
                            shear::run_2d::<CELLS>(
                                t0, tend, l, dtmin, dtmax, etaovers, zetaovers, tempcut, freezeout,
                                nb_trento,
                            );
                        }
                        ToSimulate::Trento { dt, nb_trento } => {
                            let dt = dt.unwrap_or(dx * 0.1);
                            shear::run_trento_2d::<CELLS>(
                                t0, tend, l, dt, etaovers, zetaovers, tempcut, freezeout, nb_trento,
                            );
                        }
                    }
                }
            }
        }
    }
}

fn main() {
    const STACK_SIZE: usize = 64 * 1024 * 1024; // if you want to run 2D simulation with more than 200x200 cells, you will need to increase the stack size
    thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(big_stack)
        .unwrap()
        .join()
        .unwrap();
}
