pub fn euler() -> ([[f64; 1]; 1], Option<[f64; 1]>) {
    ([[1.0]], None)
}
pub fn heun() -> ([[f64; 2]; 2], Option<[f64; 2]>) {
    ([[1.0, 0.0], [0.5, 0.5]], None)
}
pub fn radauiia2() -> ([[f64; 2]; 2], Option<[f64; 2]>) {
    ([[5.0 / 12.0, -1.0 / 12.0], [3.0 / 4.0, 1.0 / 4.0]], None)
}
pub fn gauss_legendre_1() -> ([[f64; 1]; 1], Option<[f64; 1]>) {
    ([[0.5]], Some([1.0]))
}
pub fn gauss_legendre_2() -> ([[f64; 2]; 2], Option<[f64; 2]>) {
    let sq3 = 3.0f64.sqrt();
    (
        [
            [1.0 / 4.0, 1.0 / 4.0 - 1.0 / 6.0 * sq3],
            [1.0 / 4.0 + 1.0 / 6.0 * sq3, 1.0 / 4.0],
        ],
        Some([1.0 / 2.0, 1.0 / 2.0]),
    )
}
