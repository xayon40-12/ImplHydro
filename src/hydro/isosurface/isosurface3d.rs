use std::fs::File;
use std::io::{Result, Write};

use crate::solver::context::Arr;

pub type IsoSurface3DFun<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize> =
    &'a dyn Fn(
        &Arr<C, VX, VY, VZ>, // old fields
        &Arr<C, VX, VY, VZ>, // new fields
        usize,               // ID e
        [usize; 4],          // ID [ut,ux,uy,uz]
        Option<[usize; 10]>, // ID [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
        Option<usize>,       // ID bulk
        f64,                 // time for fields
        f64,                 // dx
        f64,                 // dt
        f64,                 // freezeout energy fm^-4
    ) -> Vec<Surface3D>;

pub struct IsoSurface3DHandler<
    'a,
    const C: usize,
    const VX: usize,
    const VY: usize,
    const VZ: usize,
> {
    file: File,
    e_id: usize,
    u_ids: [usize; 4],              // [ut,ux,uy]
    shear_ids: Option<[usize; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
    bulk_id: Option<usize>,
    dx: f64,
    freezeout_energy: f64, // fm^-4
    iso_surface_fun: IsoSurface3DFun<'a, C, VX, VY, VZ>,
}

impl<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize>
    IsoSurface3DHandler<'a, C, VX, VY, VZ>
{
    pub fn new(
        filename: &str,
        e_id: usize,
        u_ids: [usize; 4],              // [ut,ux,uy,uz]
        shear_ids: Option<[usize; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
        bulk_id: Option<usize>,
        dx: f64,
        freezeout_energy: f64, // fm^-4
    ) -> Result<IsoSurface3DHandler<'a, C, VX, VY, VZ>> {
        let file = File::create(filename)?;
        let handler = IsoSurface3DHandler {
            file,
            e_id,
            u_ids,
            shear_ids,
            bulk_id,
            dx,
            freezeout_energy,
            iso_surface_fun: &zigzag3D,
        };
        Ok(handler)
    }

    pub fn find_surfaces(
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
pub struct Surface3D {
    pub pos: [f64; 4],         // [t,x,y,z]
    pub sigma: [f64; 4],       // [sigma_t, sigma_x, sigma_y, sigma_z]
    pub v: [f64; 3],           // [vx, vy, vz]
    pub pi: Option<[f64; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
    pub bulk: Option<f64>,
}

impl Surface3D {
    pub fn to_column(&self) -> [f64; 11] {
        let mut res = [0.0; 11];
        res[0..4].copy_from_slice(&self.pos);
        res[4..8].copy_from_slice(&self.sigma);
        res[8..11].copy_from_slice(&self.v);
        res
    }
    pub fn to_column_viscous(&self) -> [f64; 22] {
        let mut res = [0.0; 22];
        res[0..4].copy_from_slice(&self.pos);
        res[4..8].copy_from_slice(&self.sigma);
        res[8..11].copy_from_slice(&self.v);
        if let Some(pi) = &self.pi {
            res[11..21].copy_from_slice(pi);
        }
        if let Some(b) = &self.bulk {
            res[21] = *b;
        }
        res
    }
}

#[allow(non_snake_case)]
pub fn zigzag3D<const C: usize, const VX: usize, const VY: usize, const VZ: usize>(
    fields: &Arr<C, VX, VY, VZ>,
    new_fields: &Arr<C, VX, VY, VZ>,
    e_id: usize,
    u_ids: [usize; 4],              // [ut,ux,uy,uz]
    shear_ids: Option<[usize; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
    bulk_id: Option<usize>,
    t: f64, // time for fields
    dx: f64,
    dt: f64,               // where t+dt is the time of new_fields
    freezeout_energy: f64, // fm^-4
) -> Vec<Surface3D> {
    let mut surfaces: Vec<Surface3D> = vec![];

    let dtdx = dt * dx;
    let dxdx = dx * dx;

    let vx2 = (VX - 1) as f64 * 0.5;
    let vy2 = (VY - 1) as f64 * 0.5;
    let vz2 = (VZ - 1) as f64 * 0.5;

    for lz in 0..VZ - 1 {
        for ly in 0..VY - 1 {
            for lx in 0..VX - 1 {
                let f = fields[lz][ly][lx];
                let fx = fields[lz][ly][lx + 1];
                let fy = fields[lz][ly + 1][lx];
                let fz = fields[lz + 1][ly][lx];
                let ft = new_fields[lz][ly][lx];
                let e = f[e_id] - freezeout_energy;

                let lut = f[u_ids[0]];
                let lvx = f[u_ids[1]] / lut;
                let lvy = f[u_ids[2]] / lut;
                let lvz = f[u_ids[3]] / lut;

                let lpi = shear_ids.and_then(|ids| {
                    let mut pi = [0.0; 10];
                    for i in 0..10 {
                        pi[i] = f[ids[i]];
                    }
                    Some(pi)
                });

                let lbulk = bulk_id.and_then(|bid| Some(f[bid]));

                for (fr, ht, hx, hy, hz, s, dir) in [
                    (ft, 0.5, 0.0, 0.0, 0.0, dxdx, 0),
                    (fx, 0.0, 0.5, 0.0, 0.0, dtdx, 1),
                    (fy, 0.0, 0.0, 0.5, 0.0, dtdx, 2),
                    (fz, 0.0, 0.0, 0.0, 0.5, dtdx, 3),
                ] {
                    let er = fr[e_id] - freezeout_energy;
                    if e * er <= 0.0 {
                        let rut = fr[u_ids[0]];
                        let rvx = fr[u_ids[1]] / rut;
                        let rvy = fr[u_ids[2]] / rut;
                        let rvz = fr[u_ids[3]] / rut;

                        let t = t + ht * dt;
                        let x = (lx as f64 - vx2 + hx) * dx;
                        let y = (ly as f64 - vy2 + hy) * dx;
                        let z = (lz as f64 - vz2 + hz) * dx;

                        let vx = (lvx + rvx) * 0.5;
                        let vy = (lvy + rvy) * 0.5;
                        let vz = (lvz + rvz) * 0.5;

                        let mut sigma = [0.0; 4];
                        sigma[dir] = e.signum() * s;

                        let bulk = bulk_id.and_then(|bid| Some((fr[bid] + lbulk.unwrap()) * 0.5));
                        let pi = shear_ids.and_then(|ids| {
                            let lpi = lpi.unwrap();
                            let mut pi = [0.0; 10];
                            for i in 0..10 {
                                pi[i] = (lpi[i] + fr[ids[i]]) * 0.5;
                            }
                            Some(pi)
                        });

                        let surf = Surface3D {
                            pos: [t, x, y, z],
                            sigma,
                            v: [vx, vy, vz],
                            pi,
                            bulk,
                        };
                        surfaces.push(surf);
                    }
                }
            }
        }
    }

    surfaces
}
