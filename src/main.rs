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
    #[arg(long, default_value_t = 1.0)]
    t0: f64,
    #[arg(long, default_value_t = 4.5)]
    tend: f64,

    #[command(subcommand)]
    command: Hydro,
}

#[derive(Subcommand, Debug)]
enum Hydro {
    Ideal {
        #[command(subcommand)]
        command: ToSimulate,
    },
    Shear {
        #[arg(short, long, default_value_t = 0.08)]
        etaovers: f64,
        #[arg(long, default_value_t = 50.0)]
        tempcut: f64,
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
    let t0 = args.t0;
    let tend = args.tend;
    if t0 >= tend {
        panic!("The initial time 't0' must be smaller than the end time 'tend'.");
    }
    let checkdt = |dtmin: f64, dtmax: f64| {
        if dtmin >= dtmax {
            panic!("The smallest time interval 'dtmin={}' must be smaller than the largest 'dtmax={}'.", dtmin, dtmax);
        }
    };
    let dx = l / CELLS as f64;
    match args.command {
        Hydro::Ideal { command } => match command {
            ToSimulate::Benchmark {
                dtmin,
                dtmax,
                nb_trento,
            } => {
                let dtmax = dtmax.unwrap_or(dx * 0.1);
                checkdt(dtmin, dtmax);
                ideal::run::<CELLS>(t0, tend, l, dtmin, dtmax, nb_trento);
            }
            ToSimulate::Trento { dt, nb_trento } => {
                let dt = dt.unwrap_or(dx * 0.1);
                ideal::run_trento::<CELLS>(t0, tend, l, dt, nb_trento);
            }
        },
        Hydro::Shear {
            etaovers,
            tempcut,
            command,
        } => match command {
            ToSimulate::Benchmark {
                dtmin,
                dtmax,
                nb_trento,
            } => {
                let dtmax = dtmax.unwrap_or(dx * 0.1);
                checkdt(dtmin, dtmax);
                shear::run::<CELLS>(t0, tend, l, dtmin, dtmax, etaovers, tempcut, nb_trento);
            }
            ToSimulate::Trento { dt, nb_trento } => {
                let dt = dt.unwrap_or(dx * 0.1);
                shear::run_trento::<CELLS>(t0, tend, l, dt, etaovers, tempcut, nb_trento);
            }
        },
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
