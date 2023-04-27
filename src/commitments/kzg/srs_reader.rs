use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::{
            curves::bls12_381::{
                curve::BLS12381Curve,
                field_extension::{BLS12381PrimeField, Degree2ExtensionField},
                twist::BLS12381TwistCurve,
            },
            point::ShortWeierstrassProjectivePoint,
        },
        traits::IsEllipticCurve,
    },
    field::element::FieldElement,
};

use super::{error::KzgError, G1, G2, SRS};

macro_rules! check {
    ($cond:expr) => {
        ($cond).then_some(()).ok_or_else(|| {
            dbg!($cond);
            KzgError::BadArgs
        })
    };
    ($one:expr, $other:expr) => {
        ($one == $other).then_some(()).ok_or_else(|| {
            dbg!($one);
            dbg!($other);
            KzgError::BadArgs
        })
    };
}

// TODO: This should be an env var
const FIELD_ELEMENTS_PER_BLOB: usize = 4;
const TRUSTED_SETUP_NUM_G1_POINTS: usize = FIELD_ELEMENTS_PER_BLOB;
// TODO: change this, this only applies to Plonk
const TRUSTED_SETUP_NUM_G2_POINTS: usize = 2;
const BYTES_PER_G1: usize = 48;
const BYTES_PER_G2: usize = 96;

/// Load trusted setup from a file.
///
/// @remark The file format is `n1 n2 g1_1 g1_2 ... g1_n1 g2_1 ... g2_n2` where
///     the first two numbers are in decimal and the remainder are hexstrings
///     and any whitespace can be used as separators.
///
/// @remark See also load_trusted_setup().
/// @remark The input file will not be closed.
///
/// @param[out] out Pointer to the loaded trusted setup data
/// @param[in]  in  File handle for input
pub fn load_trusted_setup_file(input: &mut File) -> Result<SRS, KzgError> {
    let mut g1_bytes = Vec::with_capacity(TRUSTED_SETUP_NUM_G1_POINTS);
    let mut g2_bytes = Vec::with_capacity(TRUSTED_SETUP_NUM_G2_POINTS);
    let mut input_lines = BufReader::new(input).lines();

    /* Read the number of g1 points */
    let read_buf = input_lines.next().unwrap().unwrap();
    check!(
        read_buf.parse::<u64>().unwrap(),
        TRUSTED_SETUP_NUM_G1_POINTS as u64
    )?;

    /* Read the number of g2 points */
    let read_buf = input_lines.next().unwrap().unwrap();
    check!(
        read_buf.parse::<u64>().unwrap(),
        TRUSTED_SETUP_NUM_G2_POINTS as u64
    )?;

    // gi points are expressed in hex, so each byte takes 2 chars
    // Read all of the g1 points, byte by byte
    for _ in 0..g1_bytes.capacity() {
        let line = input_lines.next().unwrap().unwrap();
        check!(line.len(), BYTES_PER_G1 * 2)?;
        g1_bytes.push(line);
    }

    /* Read all of the g2 points, byte by byte */
    for _ in 0..g1_bytes.capacity() {
        let line = input_lines.next().unwrap().unwrap();
        check!(line.len(), BYTES_PER_G2 * 2)?;
        g2_bytes.push(line);
    }

    load_trusted_setup(&g1_bytes, &g2_bytes)
}

fn load_trusted_setup(g1_bytes: &[String], g2_bytes: &[String]) -> Result<SRS, KzgError> {
    let g1_points = g1_bytes
        .iter()
        .map(|g1_hex| hex_to_g1_point(g1_hex).unwrap())
        .collect::<Vec<G1>>();
    let g2_points = g2_bytes
        .iter()
        .map(|g2_hex| hex_to_g2_point(g2_hex).unwrap())
        .collect::<Vec<G2>>();
    let g2_points = [g2_points[0].clone(), g2_points[1].clone()];

    Ok(SRS::new(&g1_points, &g2_points))
}

fn hex_to_g1_point(hex: &str) -> Result<G1, KzgError> {
    let (x, y): (
        FieldElement<BLS12381PrimeField>,
        FieldElement<BLS12381PrimeField>,
    ) = uncompress_g1_point(hex)?;

    let g1_point = BLS12381Curve::create_point_from_affine(x, y)
        .ok()
        .ok_or(KzgError::BadArgs)?;

    if g1_point == ShortWeierstrassProjectivePoint::neutral_element() {
        return Ok(g1_point);
    }

    if !is_point_in_g1_subgroup(&g1_point) {
        return Err(KzgError::BadArgs);
    }

    Ok(g1_point)
}

fn uncompress_g1_point(
    _hex: &str,
) -> Result<
    (
        FieldElement<BLS12381PrimeField>,
        FieldElement<BLS12381PrimeField>,
    ),
    KzgError,
> {
    todo!()
}

fn is_point_in_g1_subgroup(_point: &G1) -> bool {
    todo!()
}

fn hex_to_g2_point(hex: &str) -> Result<G2, KzgError> {
    let (x, y): (
        FieldElement<Degree2ExtensionField>,
        FieldElement<Degree2ExtensionField>,
    ) = uncompress_g2_point(hex)?;

    BLS12381TwistCurve::create_point_from_affine(x, y)
        .ok()
        .ok_or(KzgError::BadArgs)
}

fn uncompress_g2_point(
    _hex: &str,
) -> Result<
    (
        FieldElement<Degree2ExtensionField>,
        FieldElement<Degree2ExtensionField>,
    ),
    KzgError,
> {
    todo!()
}
