use std::{
    alloc::{alloc_zeroed, Layout},
    mem::transmute,
};

pub trait Arrays<E> {}
impl<E: Copy, const N1: usize> Arrays<E> for [E; N1] {}
impl<E: Copy, const N1: usize, const N2: usize> Arrays<E> for [[E; N1]; N2] {}
impl<E: Copy, const N1: usize, const N2: usize, const N3: usize> Arrays<E> for [[[E; N1]; N2]; N3] {}
impl<E: Copy, const N1: usize, const N2: usize, const N3: usize, const N4: usize> Arrays<E>
    for [[[[E; N1]; N2]; N3]; N4]
{
}
impl<
        E: Copy,
        const N1: usize,
        const N2: usize,
        const N3: usize,
        const N4: usize,
        const N5: usize,
    > Arrays<E> for [[[[[E; N1]; N2]; N3]; N4]; N5]
{
}
impl<
        E: Copy,
        const N1: usize,
        const N2: usize,
        const N3: usize,
        const N4: usize,
        const N5: usize,
        const N6: usize,
    > Arrays<E> for [[[[[[E; N1]; N2]; N3]; N4]; N5]; N6]
{
}
impl<
        E: Copy,
        const N1: usize,
        const N2: usize,
        const N3: usize,
        const N4: usize,
        const N5: usize,
        const N6: usize,
        const N7: usize,
    > Arrays<E> for [[[[[[[E; N1]; N2]; N3]; N4]; N5]; N6]; N7]
{
}

pub fn boxarray<T: Arrays<E>, E: Copy>(e: E) -> Box<T> {
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
