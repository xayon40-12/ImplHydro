pub mod isosurface2d;
pub mod isosurface3d;

use crate::{solver::context::Arr, FLOAT};

pub enum Freezout {
    Below,
    Above,
}

pub trait IsoSurfaceHandler<const C: usize, const VX: usize, const VY: usize, const VZ: usize> {
    fn find_surfaces(
        &mut self,
        fields: &Arr<C, VX, VY, VZ>,
        new_fields: &Arr<C, VX, VY, VZ>,
        ot: FLOAT, // time for new_fields
        nt: FLOAT, // where t+dt is the time of new_fields
    ) -> Freezout;
}

pub fn toiso<'a, const C: usize, const VX: usize, const VY: usize, const VZ: usize>(
    i: impl IsoSurfaceHandler<C, VX, VY, VZ> + 'a,
) -> Box<dyn IsoSurfaceHandler<C, VX, VY, VZ> + 'a> {
    Box::new(i)
}
