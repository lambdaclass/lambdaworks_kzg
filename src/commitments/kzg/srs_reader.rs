use std::{fs::File, fmt::Display, io::Read};

use crate::{commitments::kzg::StructuredReferenceString, FIELD_ELEMENTS_PER_BLOB};

const TRUSTED_SETUP_NUM_G1_POINTS: usize = FIELD_ELEMENTS_PER_BLOB;
// TODO: change this, this only applies to Plonk
const TRUSTED_SETUP_NUM_G2_POINTS: usize = 2;
const BYTES_PER_G1: usize = 48;
const BYTES_PER_G2: usize = 96;

#[derive(Debug)]
pub enum KzgError {
    BadArgs
}

impl Display for KzgError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            KzgError::BadArgs =>
            write!(f, "The supplied data is invalid"),
        }
    }
}

///
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
///
pub fn load_trusted_setup_file<G1Point, G2Point>(input: File) -> Result<StructuredReferenceString<G1Point, G2Point>, KzgError> {
    let num_matches: usize;
    let num_matches_buf: [u8; 8];
    let g1_bytes: [u8; TRUSTED_SETUP_NUM_G1_POINTS * BYTES_PER_G1];
    let g2_bytes: [u8; TRUSTED_SETUP_NUM_G2_POINTS * BYTES_PER_G2];

    /* Read the number of g1 points */
    num_matches = input.read(&num_matches_buf);
    CHECK(num_matches == 1);
    CHECK(i == TRUSTED_SETUP_NUM_G1_POINTS);

    /* Read the number of g2 points */
    num_matches = fscanf(in, "%" SCNu64, &i);
    CHECK(num_matches == 1);
    CHECK(i == TRUSTED_SETUP_NUM_G2_POINTS);

    /* Read all of the g1 points, byte by byte */
    for (i = 0; i < TRUSTED_SETUP_NUM_G1_POINTS * BYTES_PER_G1; i++) {
        num_matches = fscanf(in, "%2hhx", &g1_bytes[i]);
        CHECK(num_matches == 1);
    }

    /* Read all of the g2 points, byte by byte */
    for (i = 0; i < TRUSTED_SETUP_NUM_G2_POINTS * BYTES_PER_G2; i++) {
        num_matches = fscanf(in, "%2hhx", &g2_bytes[i]);
        CHECK(num_matches == 1);
    }

    return load_trusted_setup(
        out,
        g1_bytes,
        TRUSTED_SETUP_NUM_G1_POINTS,
        g2_bytes,
        TRUSTED_SETUP_NUM_G2_POINTS
    );
}
