pub mod hydro;
pub mod run;
pub mod solver;

use std::{
    alloc::{alloc_zeroed, Layout},
    mem::transmute,
};

pub fn boxarray<T, E: Copy>(e: E) -> Box<T> {
    unsafe {
        let ptr = alloc_zeroed(Layout::new::<T>());
        let st = std::mem::size_of::<T>();
        let se = std::mem::size_of::<E>();
        assert!(st % se == 0);
        let n = st / se;
        let arr: *mut E = transmute(ptr);
        for i in 0..n {
            *arr.add(i) = e;
        }
        std::mem::transmute(ptr)
    }
}
