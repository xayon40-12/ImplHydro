use rayon::prelude::*;

pub mod context;
pub mod space;
pub mod time;
pub mod utils;

use {
    context::{Context, Integration},
    time::{explicit::explicit, fixpoint::fixpoint},
};

pub fn pfor2d<T: Send, const VX: usize, const VY: usize>(
    vss: &mut [[T; VX]; VY],
    f: &(dyn Fn((usize, usize, &mut T)) + Sync),
) {
    vss.par_iter_mut()
        .enumerate()
        .flat_map(|(vy, vs)| {
            vs.par_iter_mut()
                .enumerate()
                .map(move |(vx, v)| (vy, vx, v))
        })
        .for_each(f);
}

pub type Transform<'a, const F: usize, const C: usize> =
    &'a (dyn Fn(f64, [f64; F]) -> [f64; C] + Sync);

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
    schemename: &str,
    elapsed: f64,
    cost: usize,
    tsteps: usize,
    nbiter: [[usize; VX]; VY],
    integration: Integration,
) -> std::io::Result<()> {
    let t0 = context.t0;
    let tend = context.tend;
    let t = context.t;
    let dx = context.dx;
    let maxdt = context.maxdt;
    let v = context.vs;
    let dim = if VY == 1 { 1 } else { 2 };
    let foldername = &format!(
        "{}_{:?}{}d{}_{}_{}c_{:e}dt_{:e}dx",
        name, context.r.integration, dim, S, &context.r.name, VX, context.maxdt, context.dx
    );
    let mut res = format!("# t {:e}\n# cost {}\n# x y iter", t, cost);
    for f in 0..F {
        res = format!("{} {}", res, namesf[f]);
    }
    for c in 0..C {
        res = format!("{} {}", res, namesc[c]);
    }
    res = format!("{}\n", res);

    for j in 0..VY {
        for i in 0..VX {
            let y = (j as f64 - ((VY - 1) as f64) / 2.0) * dx;
            let x = (i as f64 - ((VX - 1) as f64) / 2.0) * dx;
            let mut s = format!("{:e} {:e} {}", x, y, nbiter[j][i]);
            let v = v[j][i];
            for f in 0..F {
                s = format!("{} {:e}", s, v[f]);
            }
            let vars = (context.transform)(t, (context.constraints)(t, v));
            for c in 0..C {
                s = format!("{} {:e}", s, vars[c]);
            }
            s = format!("{}\n", s);
            res = format!("{}{}", res, s);
        }
        res = format!("{}\n", res);
    }

    let dir = &format!("results/{}/{:e}", foldername, t);
    std::fs::create_dir_all(dir)?;
    std::fs::write(&format!("{}/data.txt", dir), res.as_bytes())?;
    let info = format!(
        "elapsed: {:e}\ntsteps: {}\nt0: {:e}\ntend: {:e}\nt: {:e}\ncost: {}\nnx: {}\nny: {}\ndx: {:e}\nmaxdt: {:e}\nintegration: {:?}\nscheme: {}\nname: {}\n",
        elapsed, tsteps, t0, tend, t, cost, VX, VY, dx, maxdt, integration, schemename, name,
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
    schemename: &str,
    integration: Integration,
    names: &([&str; F], [&str; C]),
) -> ([[[f64; F]; VX]; VY], f64, usize, usize) {
    // let mul = context.local_interaction[0] * context.local_interaction[1] * 2 + 1;
    let mut cost = 0.0;
    let mut tsteps = 0;
    use std::time::Instant;
    let now = Instant::now();

    let next_save = context.tend;
    let mut ctx_dt = None;
    let mut nbiter = [[1usize; VX]; VY];
    while context.t < context.tend {
        let d = next_save - context.t;
        if d < 0.0 || d > 2.0 * context.dt {
            if let Some(dt) = ctx_dt {
                context.dt = dt;
                ctx_dt = None;
            }
        } else if d > context.dt {
            ctx_dt = Some(context.dt);
            context.dt = d / 2.0;
        } else {
            ctx_dt = Some(context.dt);
            context.dt = d;
        }
        tsteps += 1;
        let (c, nbs) = match integration {
            Integration::Explicit => explicit(&mut context),
            Integration::FixPoint => fixpoint(&mut context),
        };

        cost += c;
        nbiter = nbs;
    }
    let cost = cost as usize;
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Elapsed: {:.4}", elapsed);
    let err = save(
        &context,
        names,
        name,
        schemename,
        elapsed,
        cost,
        tsteps,
        nbiter,
        integration,
    );
    let _elapsed = now.elapsed().as_secs_f64();
    eprintln!("Elapsed with save: {:.4}", _elapsed);
    match err {
        Err(e) => eprintln!("{}", e),
        Ok(()) => {}
    }
    (context.vs, context.t, cost, tsteps)
}
