use std::fs::File;
use std::io::{Result, Write};

use crate::solver::context::Arr;

use super::IsoSurfaceHandler;

pub type IsoSurface2DFun<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize> =
    &'a dyn Fn(
        &Arr<C, VX, VY, VZ>, // old fields
        &Arr<C, VX, VY, VZ>, // new fields
        usize,               // ID e
        [usize; 3],          // ID [ut,ux,uy]
        Option<[usize; 7]>,  // ID [pitt,pitx,pity,pixx,pixy,piyy,pizz]
        Option<usize>,       // ID bulk
        f64,                 // time for fields
        f64,                 // dx
        f64,                 // dt
        f64,                 // freezeout energy fm^-4
    ) -> Vec<Surface2D>;

pub struct IsoSurface2DHandler<
    'a,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> {
    file: File,
    e_id: usize,
    u_ids: [usize; 3],             // [ut,ux,uy]
    shear_ids: Option<[usize; 7]>, // [pitt,pitx,pity,pixx,pixy,piyy,pizz]
    bulk_id: Option<usize>,
    dx: f64,
    freezeout_energy: f64, // fm^-4
    iso_surface_fun: IsoSurface2DFun<'a, C, VX, VY, VZ>,
}

impl<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize>
    IsoSurface2DHandler<'a, C, VX, VY, VZ>
{
    pub fn new(
        filename: &str,
        e_id: usize,
        u_ids: [usize; 3],             // [ut,ux,uy]
        shear_ids: Option<[usize; 7]>, // [pitt,pitx,pity,pixx,pixy,piyy,pizz]
        bulk_id: Option<usize>,
        dx: f64,
        freezeout_energy: f64, // fm^-4
    ) -> Result<IsoSurface2DHandler<'a, C, VX, VY, VZ>> {
        let file = File::create(filename)?;
        let handler = IsoSurface2DHandler {
            file,
            e_id,
            u_ids,
            shear_ids,
            bulk_id,
            dx,
            freezeout_energy,
            iso_surface_fun: &zigzag2D,
        };
        Ok(handler)
    }
}

impl<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize>
    IsoSurfaceHandler<C, VX, VY, VZ> for IsoSurface2DHandler<'a, C, VX, VY, VZ>
{
    fn find_surfaces(
        &mut self,
        fields: &Arr<C, VX, VY, VZ>,
        new_fields: &Arr<C, VX, VY, VZ>,
        ot: f64, // time for new_fields
        nt: f64, // where t+dt is the time of new_fields
    ) {
        let surfaces = (self.iso_surface_fun)(
            fields,
            new_fields,
            self.e_id,
            self.u_ids,
            self.shear_ids,
            self.bulk_id,
            ot,
            self.dx,
            nt - ot,
            self.freezeout_energy,
        );
        let bytes = surfaces
            .into_iter()
            .flat_map(|ss| ss.to_column_viscous().into_iter())
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<_>>();
        if let Err(e) = self.file.write_all(&bytes) {
            eprintln!("Could not write surfaces in file: {}", e);
        }
        if let Err(e) = self.file.flush() {
            eprintln!("Could not flush file: {}", e);
        }
    }
}

#[derive(Debug)]
pub struct Surface2D {
    pub pos: [f64; 3],   // [t,x,y]
    pub sigma: [f64; 3], // [sigma_t, sigma_x, sigma_y]
    pub v: [f64; 2],     // [vx, vy]
    pub pi: Option<[f64; 7]>,
    pub bulk: Option<f64>,
}

impl Surface2D {
    pub fn to_column(&self) -> [f64; 8] {
        let mut res = [0.0; 8];
        res[0..3].copy_from_slice(&self.pos);
        res[3..6].copy_from_slice(&self.sigma);
        res[6..8].copy_from_slice(&self.v);
        res
    }
    pub fn to_column_viscous(&self) -> [f64; 16] {
        let mut res = [0.0; 16];
        res[0..3].copy_from_slice(&self.pos);
        res[3..6].copy_from_slice(&self.sigma);
        res[6..8].copy_from_slice(&self.v);
        if let Some(pi) = &self.pi {
            res[8..15].copy_from_slice(pi);
        }
        if let Some(b) = &self.bulk {
            res[15] = *b;
        }
        res
    }
}

#[allow(non_snake_case)]
pub fn zigzag2D<const C: usize, const VX: usize, const VY: usize, const VZ: usize>(
    fields: &Arr<C, VX, VY, VZ>,
    new_fields: &Arr<C, VX, VY, VZ>,
    e_id: usize,
    u_ids: [usize; 3],             // [ut,ux,uy]
    shear_ids: Option<[usize; 7]>, // [pitt,pitx,pity,pixx,pixy,piyy,pizz]
    bulk_id: Option<usize>,
    t: f64, // time for fields
    dx: f64,
    dt: f64,               // where t+dt is the time of new_fields
    freezeout_energy: f64, // fm^-4
) -> Vec<Surface2D> {
    let mut surfaces: Vec<Surface2D> = vec![];

    let dtdx = dt * dx;
    let dxdx = dx * dx;

    let vx2 = (VX - 1) as f64 * 0.5;
    let vy2 = (VY - 1) as f64 * 0.5;

    let lz: usize = 0;
    for ly in 0..VY - 1 {
        for lx in 0..VX - 1 {
            let f = fields[lz][ly][lx];
            let fx = fields[lz][ly][lx + 1];
            let fy = fields[lz][ly + 1][lx];
            let ft = new_fields[lz][ly][lx];
            let e = f[e_id] - freezeout_energy;

            let lut = f[u_ids[0]];
            let lvx = f[u_ids[1]] / lut;
            let lvy = f[u_ids[2]] / lut;

            let lpi = shear_ids.and_then(|ids| {
                let mut pi = [0.0; 7];
                for i in 0..7 {
                    pi[i] = f[ids[i]];
                }
                Some(pi)
            });

            let lbulk = bulk_id.and_then(|bid| Some(f[bid]));

            for (fr, ht, hx, hy, s, dir) in [
                (ft, 0.5, 0.0, 0.0, dxdx, 0),
                (fx, 0.0, 0.5, 0.0, dtdx, 1),
                (fy, 0.0, 0.0, 0.5, dtdx, 2),
            ] {
                let er = fr[e_id] - freezeout_energy;
                if e * er <= 0.0 {
                    let rut = fr[u_ids[0]];
                    let rvx = fr[u_ids[1]] / rut;
                    let rvy = fr[u_ids[2]] / rut;

                    let t = t + ht * dt;
                    let x = (lx as f64 - vx2 + hx) * dx;
                    let y = (ly as f64 - vy2 + hy) * dx;

                    let vx = (lvx + rvx) * 0.5;
                    let vy = (lvy + rvy) * 0.5;

                    let mut sigma = [0.0; 3];
                    sigma[dir] = e.signum() * s;

                    let bulk = bulk_id.and_then(|bid| Some((fr[bid] + lbulk.unwrap()) * 0.5));
                    let pi = shear_ids.and_then(|ids| {
                        let lpi = lpi.unwrap();
                        let mut pi = [0.0; 7];
                        for i in 0..7 {
                            pi[i] = (lpi[i] + fr[ids[i]]) * 0.5;
                        }
                        Some(pi)
                    });

                    let surf = Surface2D {
                        pos: [t, x, y],
                        sigma,
                        v: [vx, vy],
                        pi,
                        bulk,
                    };
                    surfaces.push(surf);
                }
            }
        }
    }

    surfaces
}
