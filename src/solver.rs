use crate::{
    hydro::isosurface::{
        isosurface2d::IsoSurface2DHandler, isosurface3d::IsoSurface3DHandler, toiso, Freezout,
        IsoSurfaceHandler,
    },
    solver::context::{Arr, BArr},
};
use boxarray::boxarray;
use std::collections::HashMap;

use crate::hydro::Viscosity;

use self::time::fixpoint::ErrThr;

pub mod context;
pub mod space;
pub mod time;
pub mod utils;

use {
    context::{Context, Integration},
    time::{explicit::explicit, fixpoint::fixpoint},
};

use custom_derive::custom_derive;
use enum_derive::{enum_derive_util, EnumDisplay, EnumFromStr};

// pub const EXACT: bool = false;
pub const EXACT: bool = true;
const ERROR_PROPAGATION: bool = true;
const DEBUG: bool = false;

custom_derive! {
    #[derive(Debug, Clone, Copy, EnumDisplay, EnumFromStr)]
    pub enum Solver {
        Both,
        Implicit,
        Explicit,
    }
}

pub type Constraint<'a, const F: usize, const C: usize> =
    &'a (dyn Fn(f64, [f64; F]) -> ([f64; F], [f64; C]) + Sync);
pub type Transform<'a, const F: usize> = &'a (dyn Fn(f64, [f64; F]) -> [f64; F] + Sync);
pub type Observable<
    'a,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> = (
    &'a str,
    &'a (dyn Fn(f64, &Arr<F, VX, VY, VZ>, &Arr<C, VX, VY, VZ>) -> Vec<f64> + Sync),
);

pub fn save<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const S: usize,
>(
    context: &Context<Opt, F, C, VX, VY, VZ, S>,
    (namesf, namesc): &([&str; F], [&str; C]),
    (name, case): &(&str, usize), // simulation name
    foldername: &str,
    viscosity: Viscosity,
    elapsed: f64,
    cost: usize,
    tsteps: usize,
    nbiter: &[[[usize; VX]; VY]; VZ],
    fails: usize,
    observables: &[Observable<F, C, VX, VY, VZ>],
    save_raw: bool,
) -> std::io::Result<()> {
    let t0 = context.t0;
    let tend = context.tend;
    let t = context.t;
    let dx = context.dxs[0];
    let dy = context.dxs[1];
    let dz = context.dxs[2];
    let maxdt = context.maxdt;
    let (v, trs) = &context.vstrs;
    let diffv = &context.total_diff_vs;
    let schemename = if EXACT {
        format!("{}", context.r.name)
    } else {
        format!("{}_F", context.r.name)
    };
    let integration = context.r.integration;
    let stages = S;

    let dir = &format!("{}/{:e}", foldername, t);
    std::fs::create_dir_all(dir)?;
    for (name, obs) in observables {
        let ligne = obs(t, &v, &trs);
        std::fs::write(
            &format!("{}/{}.dat", dir, name),
            ligne
                .into_iter()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        )?;
    }

    if save_raw {
        // const TXT: bool = false;
        const DIFF: bool = false;

        // let mut res = format!("# t {:e}\n# cost {}\n# x y iter", t, cost);
        // if TXT {
        //     for f in 0..F {
        //         res = format!("{} {}", res, namesf[f]);
        //     }
        //     for c in 0..C {
        //         res = format!("{} {}", res, namesc[c]);
        //     }
        //     res = format!("{}\n", res);
        // }

        let fc4 = F + C + 4;
        let mut ligne = vec![0.0f64; VY * VX * VZ * fc4];
        let f3 = F + 3;
        let mut ligne_diff = vec![0.0f64; VY * VX * VZ * f3];
        for k in 0..VZ {
            for j in 0..VY {
                for i in 0..VX {
                    let z = (k as f64 - ((VZ - 1) as f64) / 2.0) * dz;
                    let y = (j as f64 - ((VY - 1) as f64) / 2.0) * dy;
                    let x = (i as f64 - ((VX - 1) as f64) / 2.0) * dx;
                    let v = v[k][j][i];
                    let vars = trs[k][j][i];
                    let diffv = diffv[k][j][i];

                    // if TXT {
                    //     let mut s = format!("{:e} {:e} {}", x, y, nbiter[j][i]);
                    //     for f in 0..F {
                    //         s = format!("{} {:e}", s, v[f]);
                    //     }
                    //     for c in 0..C {
                    //         s = format!("{} {:e}", s, vars[c]);
                    //     }
                    //     s = format!("{}\n", s);
                    //     res = format!("{}{}", res, s);
                    // }

                    ligne[0 + fc4 * (i + VX * (j + VY * k))] = x;
                    ligne[1 + fc4 * (i + VX * (j + VY * k))] = y;
                    ligne[2 + fc4 * (i + VX * (j + VY * k))] = z;
                    ligne[3 + fc4 * (i + VX * (j + VY * k))] = nbiter[k][j][i] as f64;
                    if DIFF {
                        ligne_diff[0 + f3 * (i + VX * (j + VY * k))] = x;
                        ligne_diff[1 + f3 * (i + VX * (j + VY * k))] = y;
                        ligne_diff[2 + f3 * (i + VX * (j + VY * k))] = z;
                    }
                    for f in 0..F {
                        ligne[4 + f + fc4 * (i + VX * (j + VY * k))] = v[f];
                        if DIFF {
                            ligne_diff[3 + f + f3 * (i + VX * (j + VY * k))] = diffv[f];
                        }
                    }
                    for c in 0..C {
                        ligne[4 + F + c + fc4 * (i + VX * (j + VY * k))] = vars[c];
                    }
                }
                // if TXT {
                //     res = format!("{}\n", res);
                // }
            }
        }

        std::fs::write(
            &format!("{}/data.dat", dir),
            ligne
                .into_iter()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        )?;
        if DIFF {
            std::fs::write(
                &format!("{}/diff.dat", dir),
                ligne_diff
                    .into_iter()
                    .flat_map(|v| v.to_le_bytes())
                    .collect::<Vec<_>>(),
            )?;
        }
        // if TXT {
        //     std::fs::write(&format!("{}/data.txt", dir), res.as_bytes())?;
        // }
        let viscosity = match viscosity {
            Viscosity::Ideal => format!("Ideal"),
            Viscosity::Bulk(zeta, energycut) => {
                format!("Bulk\nzeta: {:e}\nenergycut: {:e}", zeta, energycut)
            }
            Viscosity::Shear((min, slope, crv), energycut) => {
                format!(
                    "Shear\netaovers: ({:e},{:e},{:e})\nenergycut: {:e}",
                    min, slope, crv, energycut
                )
            }
            Viscosity::Both((e_min, e_slope, e_crv), (z_max, z_width, z_peak), energycut) => {
                format!(
            "Both\netaovers: ({:e},{:e},{:e})\nzetaovers: ({:e},{:e},{:e})\nenergycut: {:e}",
            e_min, e_slope, e_crv, z_max, z_width, z_peak, energycut
        )
            }
        };
        let variables = ["x", "y", "z", "iter"]
            .iter()
            .chain(namesf.iter())
            .chain(namesc.iter())
            .map(|&s| s)
            .collect::<Vec<&str>>()
            .join(" ");
        let info = format!(
            "elapsed: {:e}\ntsteps: {}\nfails: {}\nt0: {:e}\ntend: {:e}\nt: {:e}\ncost: {}\nnx: {}\nny: {}\nnz: {}\ndx: {:e}\ndy: {:e}\ndz: {:e}\nmaxdt: {:e}\nintegration: {:?}\nscheme: {}\nstages: {}\nname: {}\ncase: {}\nviscosity: {}\nvariables: {}\n",
            elapsed, tsteps, fails, t0, tend, t, cost, VX, VY, VZ, dx, dy, dz, maxdt, integration, schemename, stages, name, case, viscosity, variables,
        );
        std::fs::write(&format!("{}/info.txt", dir), info.as_bytes())?;
    }

    Ok(())
}
pub fn save_info<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const S: usize,
>(
    context: &Context<Opt, F, C, VX, VY, VZ, S>,
    (namesf, namesc): &([&str; F], [&str; C]),
    (name, case): &(&str, usize), // simulation name
    foldername: &str,
    viscosity: Viscosity,
) -> std::io::Result<()> {
    let t0 = context.t0;
    let tend = context.tend;
    let dx = context.dxs[0];
    let dy = context.dxs[1];
    let dz = context.dxs[2];
    let maxdt = context.maxdt;
    let schemename = context.r.name;
    let integration = context.r.integration;
    let stages = S;

    let viscosity = match viscosity {
        Viscosity::Ideal => format!("Ideal"),
        Viscosity::Bulk(zeta, energycut) => {
            format!("Bulk\nzeta: {:e}\nenergycut: {:e}", zeta, energycut)
        }
        Viscosity::Shear((min, slope, crv), energycut) => {
            format!(
                "Shear\netaovers: ({:e},{:e},{:e})\nenergycut: {:e}",
                min, slope, crv, energycut
            )
        }
        Viscosity::Both((e_min, e_slope, e_crv), (z_max, z_width, z_peak), energycut) => {
            format!(
                "Both\netaovers: ({:e},{:e},{:e})\nzetaovers: ({:e},{:e},{:e})\nenergycut: {:e}",
                e_min, e_slope, e_crv, z_max, z_width, z_peak, energycut
            )
        }
    };
    let freezeout = if let Some(e_freeze) = context.freezeout_energy {
        format!("freezeoutenergy: {e_freeze:e}\n")
    } else {
        format!("")
    };
    let variables = ["x", "y", "z", "iter"]
        .iter()
        .chain(namesf.iter())
        .chain(namesc.iter())
        .map(|&s| s)
        .collect::<Vec<&str>>()
        .join(" ");
    let info = format!(
            "t0: {:e}\ntend: {:e}\nnx: {}\nny: {}\nnz: {}\ndx: {:e}\ndy: {:e}\ndz: {:e}\nmaxdt: {:e}\nintegration: {:?}\nscheme: {}\nstages: {}\nname: {}\ncase: {}\nviscosity: {}\n{}variables: {}\n",
            t0, tend, VX, VY, VZ, dx, dy, dz, maxdt, integration, schemename, stages, name, case, viscosity, freezeout, variables,
        );
    std::fs::write(&format!("{}/info.txt", foldername), info.as_bytes())?;

    Ok(())
}

pub fn run<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
    const S: usize,
>(
    mut context: Context<Opt, F, C, VX, VY, VZ, S>,
    name: &(&str, usize),
    viscosity: Viscosity,
    names: &([&str; F], [&str; C]),
    observables: &[Observable<F, C, VX, VY, VZ>],
    err_thr: ErrThr<F, C, VX, VY, VZ>,
    save_raw: Option<f64>,
) -> Option<(
    (BArr<F, VX, VY, VZ>, BArr<C, VX, VY, VZ>),
    f64,
    usize,
    usize,
)> {
    let do_save = save_raw.is_some() || observables.len() > 0;

    let dim = if VY == 1 {
        1
    } else if VZ == 1 {
        2
    } else {
        3
    };
    let schemename = if EXACT {
        format!("{}", context.r.name)
    } else {
        format!("{}_F", context.r.name)
    };
    let foldername = &format!(
        "results/{}{}_{:?}_{:?}{}d{}_{}_{}c_{:e}dt_{:e}dx{:e}dy{:e}dz",
        name.0,
        name.1,
        viscosity,
        context.r.integration,
        dim,
        S,
        &schemename,
        VX,
        context.maxdt,
        context.dxs[0],
        context.dxs[1],
        context.dxs[2],
    );
    let err = std::fs::create_dir_all(foldername);
    match err {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    let surface_filename = format!("{}/surface.dat", foldername);
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let ids = names
        .1
        .iter()
        .enumerate()
        .map(|(i, v)| (*v, i))
        .collect::<HashMap<&str, usize>>();

    let bulk_id = ids.get("Pi").and_then(|b| Some(*b));
    let mut isosurface: Option<Box<dyn IsoSurfaceHandler<C, VX, VY, VZ>>> =
        context.freezeout_energy.and_then(|freezeout_energy| {
            if VX > 1 && VY > 1 && VZ > 1 {
                let pis = [
                    "pi00", "pi01", "pi02", "pi03", "pi11", "pi12", "pi13", "pi22", "pi23", "pi33",
                ];
                let shear_ids = pis.into_iter().map(|k| ids.get(k)).enumerate().fold(
                    Some([0; 10]),
                    |acc, (i, id)| {
                        id.and_then(|id| {
                            acc.and_then(|mut pi| {
                                pi[i] = *id;
                                Some(pi)
                            })
                        })
                    },
                );
                IsoSurface3DHandler::new(
                    &surface_filename,
                    ids["e"],
                    [ids["ut"], ids["ux"], ids["uy"], ids["uz"]],
                    shear_ids,
                    bulk_id,
                    context.dxs,
                    freezeout_energy,
                )
                .map_or_else(
                    |err| {
                        eprintln!("Error in initializing IsoSurface3DHandler: {}", err);
                        None
                    },
                    |iso| Some(toiso(iso)),
                )
            } else if VX > 1 && VY > 1 {
                let pis = ["pi00", "pi01", "pi02", "pi11", "pi12", "pi22", "pi33"];
                let shear_ids = pis.into_iter().map(|k| ids.get(k)).enumerate().fold(
                    Some([0; 7]),
                    |acc, (i, id)| {
                        id.and_then(|id| {
                            acc.and_then(|mut pi| {
                                pi[i] = *id;
                                Some(pi)
                            })
                        })
                    },
                );
                IsoSurface2DHandler::new(
                    &surface_filename,
                    ids["e"],
                    [ids["ut"], ids["ux"], ids["uy"]],
                    shear_ids,
                    bulk_id,
                    context.dxs,
                    freezeout_energy,
                )
                .map_or_else(
                    |err| {
                        eprintln!("Error in initializing IsoSurface2DHandler: {}", err);
                        None
                    },
                    |iso| Some(toiso(iso)),
                )
            } else {
                None
            }
        });
    let now = Instant::now();

    let integration = context.r.integration;
    let save_every = save_raw.unwrap_or(f64::MAX).max(context.maxdt);
    let mut current_save = context.t;
    let mut next_save = current_save + save_every;
    let mut nbiter: Box<[[[usize; VX]; VY]; VZ]> = boxarray(1);
    let mut fails = 0;
    let save = |ctx: &Context<Opt, F, C, VX, VY, VZ, S>,
                cost,
                tsteps,
                nbiter: &[[[usize; VX]; VY]; VZ],
                fails| {
        let elapsed = now.elapsed().as_secs_f64();
        let err = save(
            ctx,
            names,
            name,
            foldername,
            viscosity,
            elapsed,
            cost as usize,
            tsteps,
            nbiter,
            fails,
            observables,
            save_raw.is_some(),
        );
        match err {
            Err(e) => eprintln!("{}", e),
            Ok(()) => {}
        }
    };
    match save_info(&context, names, name, &foldername, viscosity) {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    if do_save {
        save(&context, cost, tsteps, &nbiter, fails);
    }
    let m = 1e-13;
    let r = 1e10;
    let enable_save = true;
    let mut again = true;
    let mut freezout = None;
    while again {
        if do_save {
            let mut d = next_save - context.t;
            if d <= 2.0 * context.dt && enable_save {
                let mut tmp_ctx = context.clone();
                while d > tmp_ctx.t * m {
                    let n = if d <= tmp_ctx.dt { 1 } else { 2 };
                    for _ in 0..n {
                        tmp_ctx.dt = d / n as f64;
                        let res = match integration {
                            Integration::Explicit => explicit(&mut tmp_ctx),
                            Integration::FixPoint => fixpoint(&mut tmp_ctx, err_thr),
                        };
                        if res.is_none() {
                            eprintln!("Integration failed, abort current run.");
                            std::fs::remove_dir_all(foldername).expect(&format!(
                                "Could not delete directory \"{}\" after abort.",
                                foldername
                            ));
                            eprintln!("");
                            return None;
                        }
                    }
                    d = next_save - tmp_ctx.t;
                }
                tmp_ctx.t = (tmp_ctx.t * r).round() / r; // round time for saving
                save(&tmp_ctx, cost, tsteps, &nbiter, fails);
                current_save = next_save;
                next_save = (current_save + save_every).min(context.tend);
            }
        }
        tsteps += 1;
        let res = match integration {
            Integration::Explicit => explicit(&mut context),
            Integration::FixPoint => fixpoint(&mut context, err_thr),
        };
        if let Some((c, nbi, nbf)) = res {
            cost += c;
            nbiter = nbi;
            fails += nbf;
            if let Some(isosurface) = &mut isosurface {
                let (_, otrs) = &context.ovstrs;
                let (_, trs) = &context.vstrs;
                freezout = Some(isosurface.find_surfaces(otrs, trs, context.ot, context.t));
            }
        } else {
            eprintln!("Integration failed, abort current run.");
            eprintln!("");
            return None;
        }
        if context.tend <= 0.0 {
            again = freezout
                .as_ref()
                .and_then(|f| {
                    Some(match f {
                        Freezout::Below => false,
                        Freezout::Above => true,
                    })
                })
                .unwrap_or(false);
        } else {
            again = context.t < context.tend + m;
        }
    }
    let elapsed = now.elapsed().as_secs_f64();
    let hour = (elapsed / 3600.0) as u64;
    let min = ((elapsed - hour as f64 * 3600.0) / 60.0) as u64;
    let sec = (elapsed - hour as f64 * 3600.0 - min as f64 * 60.0) as u64;
    eprintln!("Elapsed: {}h{}m{}s", hour, min, sec);
    Some((context.vstrs, context.t, cost as usize, tsteps))
}
