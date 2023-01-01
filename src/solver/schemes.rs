pub type Scheme<const S: usize> = ([[f64; S]; S], Option<[f64; S]>);

pub fn euler() -> Scheme<1> {
    ([[1.0]], None)
}
pub fn heun() -> Scheme<2> {
    ([[1.0, 0.0], [0.5, 0.5]], None)
}
pub fn radauiia2() -> Scheme<2> {
    ([[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]], None)
}
pub fn gauss_legendre_1() -> Scheme<1> {
    ([[0.5]], Some([1.0]))
}
pub fn gauss_legendre_2() -> Scheme<2> {
    let sq3 = 3.0f64.sqrt();
    (
        [
            [1.0 / 4.0, 1.0 / 4.0 - 1.0 / 6.0 * sq3],
            [1.0 / 4.0 + 1.0 / 6.0 * sq3, 1.0 / 4.0],
        ],
        Some([1.0 / 2.0, 1.0 / 2.0]),
    )
}
