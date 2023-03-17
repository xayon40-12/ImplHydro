use crate::hydro::Viscosity;

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
    Opt: Sync,
    const F: usize,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const S: usize,
>(
    context: &Context<Opt, F, C, VX, VY, S>,
    (namesf, namesc): &([&str; F], [&str; C]),
    name: &str, // simulation name
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
    let dim = if VY == 1 { 1 } else { 2 };
    let schemename = context.r.name;
    let integration = context.r.integration;
    let stages = S;
    let foldername = &format!(
        "{}_{:?}_{:?}{}d{}_{}_{}c_{:e}dt_{:e}dx",
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

    let dir = &format!("results/{}/{:e}", foldername, t);
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
        Viscosity::Shear(etaovers, energycut) => {
            format!(
                "Shear\netaovers: {:e}\nenergycut: {:e}",
                etaovers, energycut
            )
        }
        Viscosity::Both(zeta, etaovers, energycut) => format!(
            "Both\nzeta: {:e}\netaovers: {:e}\nenergycut: {:e}",
            zeta, etaovers, energycut
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
    Opt: Sync,
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
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let now = Instant::now();

    let integration = context.r.integration;
    let save_every = 0.1f64.max(context.maxdt);
    let mut current_save = context.t;
    let mut next_save = current_save + save_every;
    let mut ctx_dt = None;
    let mut nbiter = [[1usize; VX]; VY];
    let mut fails = 0;
    let save = |ctx: &Context<Opt, F, C, VX, VY, S>, cost, tsteps, nbiter, fails| {
        let elapsed = now.elapsed().as_secs_f64();
        let err = save(
            ctx,
            names,
            name,
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
    while context.t < context.tend {
        let d = next_save - context.t;
        if d <= context.t * 1e-14 {
            save(&context, cost, tsteps, nbiter, fails);
            current_save = next_save;
            next_save = current_save + save_every;
            if let Some(dt) = ctx_dt {
                context.dt = dt;
                ctx_dt = None;
            }
        } else if d <= context.dt {
            ctx_dt = ctx_dt.or_else(|| Some(context.dt));
            context.dt = d;
        } else if d <= 2.0 * context.dt {
            ctx_dt = ctx_dt.or_else(|| Some(context.dt));
            context.dt = d / 2.0;
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
        } else {
            eprintln!("Integration failed, abort current run.");
            eprintln!("");
            return None;
        }
    }
    let d = next_save - context.t;
    if d < context.t * 1e-14 && context.t <= context.tend * (1.0 + 1e-14) {
        save(&context, cost, tsteps, nbiter, fails);
    }
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Elapsed: {:.4}", elapsed);
    Some((context.vstrs, context.t, cost as usize, tsteps))
}
