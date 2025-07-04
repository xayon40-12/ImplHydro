use std::fs::File;
use std::io::{Result, Write};

use crate::solver::context::{Arr, DIM};

use super::{Freezout, IsoSurfaceHandler};

pub type IsoSurface3DFun<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize> =
    &'a dyn Fn(
        &Arr<C, VX, VY, VZ>, // old fields
        &Arr<C, VX, VY, VZ>, // new fields
        usize,               // ID e
        [usize; 4],          // ID [ut,ux,uy,uz]
        Option<[usize; 10]>, // ID [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
        Option<usize>,       // ID bulk
        f64,                 // time for fields
        [f64; DIM],          // [dx, dy, dz]
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
    u_ids: [usize; 4],              // [ut,ux,uy,uz]
    shear_ids: Option<[usize; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
    bulk_id: Option<usize>,
    dxs: [f64; DIM],
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
        dxs: [f64; DIM],
        freezeout_energy: f64, // fm^-4
    ) -> Result<IsoSurface3DHandler<'a, C, VX, VY, VZ>> {
        let file = File::create(filename)?;
        let handler = IsoSurface3DHandler {
            file,
            e_id,
            u_ids,
            shear_ids,
            bulk_id,
            dxs,
            freezeout_energy,
            iso_surface_fun: &zigzag3D,
        };
        Ok(handler)
    }
}

impl<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize>
    IsoSurfaceHandler<C, VX, VY, VZ> for IsoSurface3DHandler<'a, C, VX, VY, VZ>
{
    fn find_surfaces(
        &mut self,
        fields: &Arr<C, VX, VY, VZ>,
        new_fields: &Arr<C, VX, VY, VZ>,
        ot: f64, // time for new_fields
        nt: f64, // where t+dt is the time of new_fields
    ) -> Freezout {
        let surfaces = (self.iso_surface_fun)(
            fields,
            new_fields,
            self.e_id,
            self.u_ids,
            self.shear_ids,
            self.bulk_id,
            ot,
            self.dxs,
            nt - ot,
            self.freezeout_energy,
        );
        let freezout = if surfaces.len() == 0 {
            if new_fields[0][0][0][self.e_id] > self.freezeout_energy {
                Freezout::Above
            } else {
                Freezout::Below
            }
        } else {
            Freezout::Above
        };
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

        freezout
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

/// Expect Milne coordinates
#[allow(non_snake_case)]
pub fn zigzag3D<const C: usize, const SX: usize, const SY: usize, const SETA: usize>(
    fields: &Arr<C, SX, SY, SETA>,
    new_fields: &Arr<C, SX, SY, SETA>,
    e_id: usize,
    u_ids: [usize; 4],              // [utau,ux,uy,ueta]
    shear_ids: Option<[usize; 10]>, // [pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz]
    bulk_id: Option<usize>,
    tau: f64, // time for fields
    dxs: [f64; DIM],
    dtau: f64,             // where t+dt is the time of new_fields
    freezeout_energy: f64, // fm^-4
) -> Vec<Surface3D> {
    let mut surfaces: Vec<Surface3D> = vec![];

    let dx = dxs[0];
    let dy = dxs[1];
    let deta = tau * dxs[2];

    let vx2 = (SX - 1) as f64 * 0.5;
    let vy2 = (SY - 1) as f64 * 0.5;
    let veta2 = (SETA - 1) as f64 * 0.5;

    for leta in 0..SETA - 1 {
        for ly in 0..SY - 1 {
            for lx in 0..SX - 1 {
                let f = fields[leta][ly][lx];
                let fx = fields[leta][ly][lx + 1];
                let fy = fields[leta][ly + 1][lx];
                let feta = fields[leta + 1][ly][lx];
                let ftau = new_fields[leta][ly][lx];
                let e = f[e_id] - freezeout_energy;

                let lutau = f[u_ids[0]];
                let lvx = f[u_ids[1]] / lutau;
                let lvy = f[u_ids[2]] / lutau;
                let lveta = f[u_ids[3]] / lutau;

                let lpi = shear_ids.and_then(|ids| {
                    let mut pi = [0.0; 10];
                    for i in 0..10 {
                        pi[i] = f[ids[i]];
                    }
                    Some(pi)
                });

                let lbulk = bulk_id.and_then(|bid| Some(f[bid]));

                // output in covariant form, so '-' for the space components
                for (fr, htau, hx, hy, heta, s, dir) in [
                    (ftau, 0.5, 0.0, 0.0, 0.0, dx * dy * deta, 0),
                    (fx, 0.0, 0.5, 0.0, 0.0, -dtau * dy * deta, 1),
                    (fy, 0.0, 0.0, 0.5, 0.0, -dtau * dx * deta, 2),
                    (feta, 0.0, 0.0, 0.0, 0.5, -dtau * dx * dy, 3),
                ] {
                    let er = fr[e_id] - freezeout_energy;
                    if e * er <= 0.0 {
                        let rutau = fr[u_ids[0]];
                        let rvx = fr[u_ids[1]] / rutau;
                        let rvy = fr[u_ids[2]] / rutau;
                        let rveta = fr[u_ids[3]] / rutau;

                        let tau = tau + htau * dtau;
                        let x = (lx as f64 - vx2 + hx) * dx;
                        let y = (ly as f64 - vy2 + hy) * dx;
                        let eta = (leta as f64 - veta2 + heta) * dx;

                        let vx = (lvx + rvx) * 0.5;
                        let vy = (lvy + rvy) * 0.5;
                        let veta = (lveta + rveta) * 0.5;

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
                            pos: [tau, x, y, eta],
                            sigma,
                            v: [vx, vy, veta],
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
