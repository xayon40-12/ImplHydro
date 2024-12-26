use clap::{command, Parser, Subcommand};
use impl_hydro::{
    run::{ideal1d, ideal2d, ideal3d, viscous2d, viscous3d},
    solver::Solver,
    FLOAT,
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(short, long, default_value_t = 100)]
    cells: usize,

    #[arg(short, long, default_value_t = 40.0)]
    physical_length: FLOAT,

    #[arg(short, long, default_value_t = 5020.0)] // GeV
    energy: FLOAT,

    #[arg(short, long, default_value_t = Solver::Both)]
    solver_type: Solver,

    #[arg(short, long)]
    raw_data_every: Option<FLOAT>,

    #[command(subcommand)]
    command: Dim,
}

#[derive(Subcommand, Debug)]
enum Dim {
    Dim1 {
        #[arg(long, default_value_t = 0.0)]
        t0: FLOAT,
        #[arg(long, default_value_t = 15.0)]
        tend: FLOAT,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim2 {
        #[arg(long, default_value_t = 1.0)]
        t0: FLOAT,
        #[arg(long, default_value_t = 0.0)]
        tend: FLOAT,
        #[command(subcommand)]
        command: Hydro,
    },
    Dim3 {
        #[arg(long, default_value_t = 1.0)]
        t0: FLOAT,
        #[arg(long, default_value_t = 0.0)]
        tend: FLOAT,
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
        etaovers_min: FLOAT,
        #[arg(long, default_value_t = 1.6)] // GeV^-1
        etaovers_slope: FLOAT,
        #[arg(long, default_value_t = -0.29)]
        etaovers_crv: FLOAT,
        #[arg(long, default_value_t = 0.032)]
        zetaovers_max: FLOAT,
        #[arg(long, default_value_t = 0.024)] // GeV
        zetaovers_width: FLOAT,
        #[arg(long, default_value_t = 0.175)] // GeV
        zetaovers_peak: FLOAT,
        #[arg(long, default_value_t = 0.020)] // GeV
        tempcut: FLOAT,
        #[arg(long, default_value_t = 0.151)] // GeV
        freezeout: FLOAT,
        #[command(subcommand)]
        command: ToSimulate,
    },
}

#[derive(Subcommand, Debug)]
enum ToSimulate {
    Benchmark {
        #[arg(long, default_value_t = 1e-3)]
        dtmin: FLOAT,
        #[arg(long)]
        dtmax: Option<FLOAT>,
        #[arg(short, long, default_value_t = 0.1)]
        dt_over_dx_max: FLOAT,
        #[arg(short, long, default_value_t = 1)]
        nb_trento: usize,
        #[arg(short, long, default_value_t = 0)]
        first_trento: usize,
    },
    Trento {
        #[arg(long)]
        dt: Option<FLOAT>,
        #[arg(short, long, default_value_t = 0.1)]
        dt_over_dx: FLOAT,
        #[arg(short, long, default_value_t = 1)]
        nb_trento: usize,
        #[arg(short, long, default_value_t = 0)]
        first_trento: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.cells {
        100 => run::<100, 25>(cli),
        200 => run::<200, 51>(cli),
        _ => panic!("The number of cells must be a value from the list {{100,200}}."),
    };
}
fn run<const XY: usize, const Z: usize>(args: Cli) {
    let xy_len = args.physical_length;
    // let sqrts = args.collision_energy;
    let etas_len = xy_len / 2.0; //2.0 * (0.5 * sqrts / 0.2).acosh();
    let checkdt = |dtmin: FLOAT, dtmax: FLOAT| {
        if dtmin >= dtmax {
            panic!("The smallest time interval 'dtmin={}' must be smaller than the largest 'dtmax={}'.", dtmin, dtmax);
        }
    };
    let checkt = |_t0: FLOAT, _tend: FLOAT| {
        // if t0 > tend {
        //     panic!("The initial time 't0' must be smaller than the end time 'tend'.");
        // }
    };
    let dx = xy_len / XY as FLOAT;
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
                        dt_over_dx_max,
                        nb_trento: _,
                        first_trento: _,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dt_over_dx_max);
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
                        dt_over_dx_max,
                        nb_trento,
                        first_trento,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dt_over_dx_max);
                        checkdt(dtmin, dtmax);
                        let nf_trento = (nb_trento, first_trento);
                        ideal2d::run_2d::<XY>(
                            solver,
                            t0,
                            tend,
                            xy_len,
                            dtmin,
                            dtmax,
                            enable_gubser,
                            nf_trento,
                            save_raw,
                        );
                    }
                    ToSimulate::Trento {
                        dt,
                        dt_over_dx,
                        nb_trento,
                        first_trento,
                    } => {
                        let dt = dt.unwrap_or(dx * dt_over_dx);
                        let nf_trento = (nb_trento, first_trento);
                        ideal2d::run_trento_2d::<XY>(
                            solver, t0, tend, xy_len, dt, nf_trento, save_raw,
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
                            dt_over_dx_max,
                            nb_trento,
                            first_trento,
                        } => {
                            let dtmax = dtmax.unwrap_or(dx * dt_over_dx_max);
                            checkdt(dtmin, dtmax);
                            let nf_trento = (nb_trento, first_trento);
                            viscous2d::run_2d::<XY>(
                                solver, t0, tend, xy_len, dtmin, dtmax, etaovers, zetaovers,
                                tempcut, freezeout, nf_trento, save_raw,
                            );
                        }
                        ToSimulate::Trento {
                            dt,
                            dt_over_dx,
                            nb_trento,
                            first_trento,
                        } => {
                            let dt = dt.unwrap_or(dx * dt_over_dx);
                            let nf_trento = (nb_trento, first_trento);
                            viscous2d::run_trento_2d::<XY>(
                                solver, t0, tend, xy_len, dt, etaovers, zetaovers, tempcut,
                                freezeout, nf_trento, save_raw,
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
                        dt_over_dx_max,
                        nb_trento,
                        first_trento,
                    } => {
                        let dtmax = dtmax.unwrap_or(dx * dt_over_dx_max);
                        checkdt(dtmin, dtmax);
                        let nf_trento = (nb_trento, first_trento);
                        ideal3d::run_3d::<XY, Z>(
                            solver, t0, tend, xy_len, etas_len, dtmin, dtmax, nf_trento, save_raw,
                        );
                    }
                    ToSimulate::Trento {
                        dt,
                        dt_over_dx,
                        nb_trento,
                        first_trento,
                    } => {
                        let dt = dt.unwrap_or(dx * dt_over_dx);
                        let nf_trento = (nb_trento, first_trento);
                        ideal3d::run_trento_3d::<XY, Z>(
                            solver, t0, tend, xy_len, etas_len, dt, nf_trento, save_raw,
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
                            dt_over_dx_max,
                            nb_trento,
                            first_trento,
                        } => {
                            let dtmax = dtmax.unwrap_or(dx * dt_over_dx_max);
                            checkdt(dtmin, dtmax);
                            let nf_trento = (nb_trento, first_trento);
                            viscous3d::run_3d::<XY, Z>(
                                solver, t0, tend, xy_len, etas_len, dtmin, dtmax, etaovers,
                                zetaovers, tempcut, freezeout, nf_trento, save_raw,
                            );
                        }
                        ToSimulate::Trento {
                            dt,
                            dt_over_dx,
                            nb_trento,
                            first_trento,
                        } => {
                            let dt = dt.unwrap_or(dx * dt_over_dx);
                            let nf_trento = (nb_trento, first_trento);
                            viscous3d::run_trento_3d::<XY, Z>(
                                solver, t0, tend, xy_len, etas_len, dt, etaovers, zetaovers,
                                tempcut, freezeout, nf_trento, save_raw,
                            );
                        }
                    }
                }
            }
        }
    }
}
