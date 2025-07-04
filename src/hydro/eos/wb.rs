// Coefficients taken from CUDA code https://github.com/bazow/gpu-vh/blob/master/rhic/rhic-trunk/src/main/cuda/edu/osu/rhic/trunk/eos/EquationOfState.cu

// Takes energy [fm^-4] and return speed of sound squared
pub fn dpde(e: f64) -> f64 {
    let a_arr = [
        5.191934309650155e-32,
        4.123605749683891e-23,
        3.1955868410879504e-16,
        1.4170364808063119e-10,
        6.087136671592452e-6,
        0.02969737949090831,
        15.382615282179595,
        460.6487249985994,
        1612.4245252438795,
        275.0492627924299,
        58.60283714484669,
        6.504847576502024,
        0.03009027913262399,
        8.189430244031285e-6,
    ];
    let b_arr = [
        1.4637868900982493e-30,
        6.716598285341542e-22,
        3.5477700458515908e-15,
        1.1225580509306008e-9,
        0.00003551782901018317,
        0.13653226327408863,
        60.85769171450653,
        1800.5461219450308,
        15190.225535036281,
        590.2572000057821,
        293.99144775704605,
        21.461303090563028,
        0.09301685073435291,
        0.000024810902623582917,
    ];

    from_arr(e, a_arr, b_arr)
}

#[test]
pub fn test_wb() {
    let n = 10000;
    let m = 1.0;
    for i in 0..n {
        let e = i as f64 / n as f64 * m;
        let t = T(e);
        let h = p(e);
        let d = dpde(e);
        println!("@ {:e} {:e} {:e} {:e}", e, t, h, d);
    }
}

// Takes energy [fm^-4] and return temperature in [fm^-1]
#[allow(non_snake_case)]
pub fn T(e: f64) -> f64 {
    let a_arr = [
        1.510073201405604e-29,
        8.014062800678687e-18,
        2.4954778310451065e-10,
        0.000063810382643387,
        0.4873490574161924,
        207.48582344326206,
        6686.07424325115,
        14109.766109389702,
        1471.6180520527757,
        14.055788949565482,
        0.015421252394182246,
        1.5780479034557783e-6,
    ];
    let b_arr = [
        7.558667139355393e-28,
        1.3686372302041508e-16,
        2.998130743142826e-9,
        0.0005036835870305458,
        2.316902328874072,
        578.0778724946719,
        11179.193315394154,
        17965.67607192861,
        1051.0730543534657,
        5.916312075925817,
        0.003778342768228011,
        1.8472801679382593e-7,
    ];

    from_arr(e, a_arr, b_arr)
}

// Takes energy [fm^-4] and return pressure in [fm^-4]
pub fn p(e: f64) -> f64 {
    let a_arr = [
        -0.25181736420168666,
        9737.845799644809,
        1.077580993288114e6,
        3.1729694865420084e6,
        1.6357487344679043e6,
        334334.4309240126,
        41913.439282708554,
        6340.448389300905,
        141.5073484468774,
        0.7158279081255019,
        0.0009417586777847889,
        3.1188455176941583e-7,
        1.9531729608963267e-11,
    ];
    let b_arr = [
        45829.44617893836,
        4.0574329080826794e6,
        2.0931169138134286e7,
        1.3512402226067686e7,
        1.7851642641834426e6,
        278581.2989342773,
        26452.34905933697,
        499.04919730607065,
        2.3405487982094204,
        0.002962497695527404,
        9.601103399348206e-7,
        5.928138360995685e-11,
        3.2581066229887368e-18,
    ];

    from_arr(e, a_arr, b_arr) + 0.25182 / 4.5829e4 // compensate negative values
}

pub fn s(e: f64) -> f64 {
    (e + p(e)) / T(e)
}

pub fn from_arr<const N: usize>(e: f64, a_arr: [f64; N], b_arr: [f64; N]) -> f64 {
    let mut a = a_arr[0];
    let mut b = b_arr[0];
    let mut ei = 1.0;
    for i in 1..N {
        ei *= e;
        a += a_arr[i] * ei;
        b += b_arr[i] * ei;
    }
    a / b
}
