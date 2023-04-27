use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::commitments::kzg::StructuredReferenceString;

use super::error::KzgError;

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
pub fn load_trusted_setup_file<G1Point, G2Point>(
    input: &mut File,
) -> Result<StructuredReferenceString<G1Point, G2Point>, KzgError> {
    let mut g1_bytes = Vec::with_capacity(TRUSTED_SETUP_NUM_G1_POINTS);
    let mut g2_bytes = Vec::with_capacity(TRUSTED_SETUP_NUM_G2_POINTS);
    let mut input_lines = BufReader::new(input).lines();

    /* Read the number of g1 points */
    let read_buf = input_lines.next().unwrap().unwrap();
    check!(
        u64::from_str_radix(&read_buf, 10).unwrap(),
        TRUSTED_SETUP_NUM_G1_POINTS as u64
    )?;

    /* Read the number of g2 points */
    let read_buf = input_lines.next().unwrap().unwrap();
    check!(
        u64::from_str_radix(&read_buf, 10).unwrap(),
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

fn load_trusted_setup<G1Point, G2Point>(
    _g1_bytes: &[String],
    _g2_bytes: &[String],
) -> Result<StructuredReferenceString<G1Point, G2Point>, KzgError> {
    todo!()
}
