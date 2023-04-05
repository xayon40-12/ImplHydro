use std::collections::HashMap;

use crate::hydro::{isosurface::IsoSurfaceHandler, Viscosity};

pub mod context;
pub mod space;
pub mod time;
pub mod utils;

use {
    context::{Context, Integration},
    time::{explicit::explicit, fixpoint::fixpoint},
};

pub type Constraint<'a, const F: usize, const C: usize> =
    &'a (dyn Fn(f64, [f64; F], [f64; C]) -> ([f64; F], [f64; C]) + Sync); // takes the new [f64; C] and the old [f64; C]
pub type Transform<'a, const F: usize> = &'a (dyn Fn(f64, [f64; F]) -> [f64; F] + Sync);
pub type Observable<'a, const F: usize, const C: usize, const VX: usize, const VY: usize> = (
    &'a str,
    &'a (dyn Fn(f64, &[[[f64; F]; VX]; VY], &[[[f64; C]; VX]; VY]) -> Vec<f64> + Sync),
);

pub fn save<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    context: &Context<Opt, F, C, VX, VY, S>,
    (namesf, namesc): &([&str; F], [&str; C]),
    name: &str, // simulation name
    foldername: &str,
    viscosity: Viscosity,
    elapsed: f64,
    cost: usize,
    tsteps: usize,
    nbiter: [[usize; VX]; VY],
    fails: usize,
    observables: &[Observable<F, C, VX, VY>],
) -> std::io::Result<()> {
    let t0 = context.t0;
    let tend = context.tend;
    let t = context.t;
    let dx = context.dx;
    let maxdt = context.maxdt;
    let (v, trs) = context.vstrs;
    let diffv = context.total_diff_vs;
    let schemename = context.r.name;
    let integration = context.r.integration;
    let stages = S;

    const TXT: bool = false;
    const DIFF: bool = false;

    let mut res = format!("# t {:e}\n# cost {}\n# x y iter", t, cost);
    if TXT {
        for f in 0..F {
            res = format!("{} {}", res, namesf[f]);
        }
        for c in 0..C {
            res = format!("{} {}", res, namesc[c]);
        }
        res = format!("{}\n", res);
    }

    let fc3 = F + C + 3;
    let mut ligne = vec![0.0f64; VY * VX * fc3];
    let f2 = F + 2;
    let mut ligne_diff = vec![0.0f64; VY * VX * f2];
    for j in 0..VY {
        for i in 0..VX {
            let y = (j as f64 - ((VY - 1) as f64) / 2.0) * dx;
            let x = (i as f64 - ((VX - 1) as f64) / 2.0) * dx;
            let v = v[j][i];
            let vars = trs[j][i];
            let diffv = diffv[j][i];

            if TXT {
                let mut s = format!("{:e} {:e} {}", x, y, nbiter[j][i]);
                for f in 0..F {
                    s = format!("{} {:e}", s, v[f]);
                }
                for c in 0..C {
                    s = format!("{} {:e}", s, vars[c]);
                }
                s = format!("{}\n", s);
                res = format!("{}{}", res, s);
            }

            ligne[0 + fc3 * (i + j * VX)] = x;
            ligne[1 + fc3 * (i + j * VX)] = y;
            ligne[2 + fc3 * (i + j * VX)] = nbiter[j][i] as f64;
            if DIFF {
                ligne_diff[0 + f2 * (i + j * VX)] = x;
                ligne_diff[1 + f2 * (i + j * VX)] = y;
            }
            for f in 0..F {
                ligne[3 + f + fc3 * (i + j * VX)] = v[f];
                if DIFF {
                    ligne_diff[2 + f + f2 * (i + j * VX)] = diffv[f];
                }
            }
            for c in 0..C {
                ligne[3 + F + c + fc3 * (i + j * VX)] = vars[c];
            }
        }
        if TXT {
            res = format!("{}\n", res);
        }
    }

    let dir = &format!("{}/{:e}", foldername, t);
    std::fs::create_dir_all(dir)?;
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
    if TXT {
        std::fs::write(&format!("{}/data.txt", dir), res.as_bytes())?;
    }
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
        Viscosity::Both((e_min, e_slope, e_crv), (z_max, z_width, z_peak), energycut) => format!(
            "Both\netaovers: ({:e},{:e},{:e})\nzetaovers: ({:e},{:e},{:e})\nenergycut: {:e}",
            e_min, e_slope, e_crv, z_max, z_width, z_peak, energycut
        ),
    };
    let variables = ["x", "y", "iter"]
        .iter()
        .chain(namesf.iter())
        .chain(namesc.iter())
        .map(|&s| s)
        .collect::<Vec<&str>>()
        .join(" ");
    let info = format!(
        "elapsed: {:e}\ntsteps: {}\nfails: {}\nt0: {:e}\ntend: {:e}\nt: {:e}\ncost: {}\nnx: {}\nny: {}\ndx: {:e}\nmaxdt: {:e}\nintegration: {:?}\nscheme: {}\nstages: {}\nname: {}\nviscosity: {}\nvariables: {}\n",
        elapsed, tsteps, fails, t0, tend, t, cost, VX, VY, dx, maxdt, integration, schemename, stages, name, viscosity, variables,
    );
    std::fs::write(&format!("{}/info.txt", dir), info.as_bytes())?;

    Ok(())
}

pub fn run<
    Opt: Clone + Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    mut context: Context<Opt, F, C, VX, VY, S>,
    name: &str,
    viscosity: Viscosity,
    names: &([&str; F], [&str; C]),
    observables: &[Observable<F, C, VX, VY>],
) -> Option<(
    ([[[f64; F]; VX]; VY], [[[f64; C]; VX]; VY]),
    f64,
    usize,
    usize,
)> {
    let dim = if VY == 1 { 1 } else { 2 };
    let foldername = &format!(
        "results/{}_{:?}_{:?}{}d{}_{}_{}c_{:e}dt_{:e}dx",
        name,
        viscosity,
        context.r.integration,
        dim,
        S,
        &context.r.name,
        VX,
        context.maxdt,
        context.dx
    );
    let err = std::fs::create_dir_all(foldername);
    match err {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    let surface_filename = format!("{}/surface.dat", foldername);
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let ids = names
        .1
        .iter()
        .enumerate()
        .map(|(i, v)| (*v, i))
        .collect::<HashMap<&str, usize>>();

    let pis = ["pi00", "pi01", "pi02", "pi11", "pi12", "pi22", "pi33"];
    let shear_ids =
        pis.into_iter()
            .map(|k| ids.get(k))
            .enumerate()
            .fold(Some([0; 7]), |acc, (i, id)| {
                id.and_then(|id| {
                    acc.and_then(|mut pi| {
                        pi[i] = *id;
                        Some(pi)
                    })
                })
            });
    let bulk_id = ids.get("Pi").and_then(|b| Some(*b));
    let mut isosurface: Option<IsoSurfaceHandler<C, VX, VY>> =
        context.freezeout_energy.and_then(|freezeout_energy| {
            IsoSurfaceHandler::new(
                &surface_filename,
                ids["e"],
                [ids["ut"], ids["ux"], ids["uy"]],
                shear_ids,
                bulk_id,
                context.dx,
                freezeout_energy,
            )
            .map_or_else(
                |err| {
                    eprintln!("Error in initializing IsoSurfaceHondler: {}", err);
                    None
                },
                |iso| Some(iso),
            )
        });
    let now = Instant::now();

    let integration = context.r.integration;
    let save_every = 0.1f64; //.max(context.maxdt);
    let mut current_save = context.t;
    let mut next_save = current_save + save_every;
    let mut nbiter = [[1usize; VX]; VY];
    let mut fails = 0;
    let save = |ctx: &Context<Opt, F, C, VX, VY, S>, cost, tsteps, nbiter, fails| {
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
        );
        match err {
            Err(e) => eprintln!("{}", e),
            Ok(()) => {}
        }
    };
    save(&context, cost, tsteps, nbiter, fails);
    let m = 1e-13;
    let r = 1e10;
    while context.t < context.tend - m {
        let mut d = next_save.min(context.tend) - context.t;
        context.dt = context.dt.min(context.tend - context.t);
        if d <= 2.0 * context.dt {
            let mut tmp_ctx = context.clone();
            while d > tmp_ctx.t * m {
                let n = if d <= tmp_ctx.dt { 1 } else { 2 };
                for _ in 0..n {
                    tmp_ctx.dt = d / n as f64;
                    match integration {
                        Integration::Explicit => explicit(&mut tmp_ctx),
                        Integration::FixPoint(err_ref_p) => fixpoint(&mut tmp_ctx, err_ref_p),
                    };
                }
                d = next_save.min(tmp_ctx.tend) - tmp_ctx.t;
            }
            tmp_ctx.t = (tmp_ctx.t * r).round() / r; // round time for saving
            save(&tmp_ctx, cost, tsteps, nbiter, fails);
            current_save = next_save;
            next_save = current_save + save_every;
        }
        tsteps += 1;
        let res = match integration {
            Integration::Explicit => explicit(&mut context),
            Integration::FixPoint(err_ref_p) => fixpoint(&mut context, err_ref_p),
        };
        if let Some((c, nbi, nbf)) = res {
            cost += c;
            nbiter = nbi;
            fails += nbf;
            if let Some(isosurface) = &mut isosurface {
                let (_, otrs) = &context.ovstrs;
                let (_, trs) = &context.vstrs;
                isosurface.find_surfaces(otrs, trs, context.ot, context.t);
            }
        } else {
            eprintln!("Integration failed, abort current run.");
            eprintln!("");
            return None;
        }
    }
    let d = next_save.min(context.tend) - context.t;
    if d < context.t * m && context.t <= context.tend * (1.0 + m) {
        context.t = (context.t * r).round() / r; // round time for saving
        save(&context, cost, tsteps, nbiter, fails);
    }
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Elapsed: {:.4}", elapsed);
    Some((context.vstrs, context.t, cost as usize, tsteps))
}
