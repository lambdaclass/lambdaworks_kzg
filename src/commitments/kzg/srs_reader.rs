use std::{fs::File, io::Read};

use crate::commitments::kzg::StructuredReferenceString;

use super::error::KzgError;

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
    let mut num_bytes_read: usize;
    let mut read_buf = [0_u8; 8];
    let mut gi_byte_hex = [0_u8; 2];
    let mut g1_bytes = [0_u8; TRUSTED_SETUP_NUM_G1_POINTS * BYTES_PER_G1];
    let mut g2_bytes = [0_u8; TRUSTED_SETUP_NUM_G2_POINTS * BYTES_PER_G2];

    /* Read the number of g1 points */
    num_bytes_read = read_file(input, &mut read_buf)?;
    check(
        num_bytes_read == 1 && u64::from_be_bytes(read_buf) == TRUSTED_SETUP_NUM_G1_POINTS as u64,
    )?;

    /* Read the number of g2 points */
    num_bytes_read = read_file(input, &mut read_buf)?;
    check(
        num_bytes_read == 1 && u64::from_be_bytes(read_buf) == TRUSTED_SETUP_NUM_G2_POINTS as u64,
    )?;

    /* Read all of the g1 points, byte by byte */
    for g1_byte in g1_bytes.iter_mut() {
        num_bytes_read = read_file(input, &mut gi_byte_hex)?;
        check(num_bytes_read == 2)?;
        *g1_byte = hex_bytes_to_u8(gi_byte_hex)?;
    }

    /* Read all of the g2 points, byte by byte */
    for g2_byte in g2_bytes.iter_mut() {
        num_bytes_read = read_file(input, &mut gi_byte_hex)?;
        check(num_bytes_read == 2)?;
        *g2_byte = hex_bytes_to_u8(gi_byte_hex)?;
    }

    load_trusted_setup(&g1_bytes, &g2_bytes)
}

fn load_trusted_setup<G1Point, G2Point>(
    _g1_bytes: &[u8],
    _g2_bytes: &[u8],
) -> Result<StructuredReferenceString<G1Point, G2Point>, KzgError> {
    todo!()
}

fn read_file(file: &mut File, buf: &mut [u8]) -> Result<usize, KzgError> {
    match file.read(buf) {
        Ok(num_bytes) => Ok(num_bytes),
        Err(_) => Err(KzgError::BadArgs),
    }
}

fn check(condition: bool) -> Result<(), KzgError> {
    condition.then_some(()).ok_or(KzgError::BadArgs)
}

fn hex_bytes_to_u8(hex: [u8; 2]) -> Result<u8, KzgError> {
    String::from_utf8(hex.to_vec())
        .ok()
        .and_then(|hex_str| u8::from_str_radix(&hex_str, 16).ok())
        .ok_or(KzgError::BadArgs)
}
