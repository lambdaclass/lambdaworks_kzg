use crate::KZGSettings;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines};
use std::path::Path;

/// Helper function that reads a file line by line and
/// returns an iterator over the lines.
fn read_lines<P>(path: P) -> io::Result<Lines<BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(path)?;
    Ok(BufReader::new(file).lines())
}

/// Load trusted setup from a file.
///
/// The file format is `n1 n2 g1_1 g1_2 ... g1_n1 g2_1 ... g2_n2` where
/// the first two numbers are in decimal and the remainder are hexstrings
/// and any whitespace can be used as separators.
///
/// See also `load_trusted_setup`.
///
/// # Arguments
///
/// * `path` - Path to the file containing the trusted setup
///
/// # Returns
///
/// * `KZGSettings` - The loaded trusted setup data
pub fn load_trusted_setup_file(path: &str) -> io::Result<u8> {
    let mut lines = read_lines(path)?;

    /*
    int num_matches;
    uint64_t i;
    uint8_t g1_bytes[TRUSTED_SETUP_NUM_G1_POINTS * BYTES_PER_G1];
    uint8_t g2_bytes[TRUSTED_SETUP_NUM_G2_POINTS * BYTES_PER_G2];
    */

    // Read the number of g1 points
    let num_g1_points = lines
        .next()
        .ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format",
        ))??
        .parse::<usize>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    let num_g2_points = lines
        .next()
        .ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format",
        ))??
        .parse::<u64>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    println!("num_g1_points: {num_g1_points}, num_g2_points: {num_g2_points}");

    for line in lines {
        if let Ok(line) = line {
            println!("{}", line);
        }
    }

    /*
    // Read all of the g1 points, byte by byte
    for (i = 0; i < TRUSTED_SETUP_NUM_G1_POINTS * BYTES_PER_G1; i++) {
        num_matches = fscanf(in, "%2hhx", &g1_bytes[i]);
        CHECK(num_matches == 1);
    }

    // Read all of the g2 points, byte by byte
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
    */

    //todo!("Implement load_trusted_setup_file")
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::load_trusted_setup_file;

    #[test]
    fn test_read_srs() {
        load_trusted_setup_file("test/trusted_setup_4.txt").unwrap();
    }
}
/*
C_KZG_RET load_trusted_setup(
    KZGSettings *out,
    const uint8_t *g1_bytes, /* n1 * 48 bytes */
    size_t n1,
    const uint8_t *g2_bytes, /* n2 * 96 bytes */
    size_t n2
);


/**
 * Load trusted setup into a KZGSettings.
 *
 * @remark Free after use with free_trusted_setup().
 *
 * @param[out] out      Pointer to the stored trusted setup data
 * @param[in]  g1_bytes Array of G1 points
 * @param[in]  n1       Number of `g1` points in g1_bytes
 * @param[in]  g2_bytes Array of G2 points
 * @param[in]  n2       Number of `g2` points in g2_bytes
 */
C_KZG_RET load_trusted_setup(
    KZGSettings *out,
    const uint8_t *g1_bytes,
    size_t n1,
    const uint8_t *g2_bytes,
    size_t n2
) {
    uint64_t i;
    blst_p2_affine g2_affine;
    g1_t *g1_projective = NULL;
    C_KZG_RET ret;

    out->fs = NULL;
    out->g1_values = NULL;
    out->g2_values = NULL;

    /* Sanity check in case this is called directly */
    CHECK(n1 == TRUSTED_SETUP_NUM_G1_POINTS);
    CHECK(n2 == TRUSTED_SETUP_NUM_G2_POINTS);

    /* Allocate all of our arrays */
    ret = new_g1_array(&out->g1_values, n1);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_g2_array(&out->g2_values, n2);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_g1_array(&g1_projective, n1);
    if (ret != C_KZG_OK) goto out_error;

    /* Convert all g1 bytes to g1 points */
    for (i = 0; i < n1; i++) {
        ret = validate_kzg_g1(
            &g1_projective[i], (Bytes48 *)&g1_bytes[BYTES_PER_G1 * i]
        );
        if (ret != C_KZG_OK) goto out_error;
    }

    /* Convert all g2 bytes to g2 points */
    for (i = 0; i < n2; i++) {
        blst_p2_uncompress(&g2_affine, &g2_bytes[BYTES_PER_G2 * i]);
        blst_p2_from_affine(&out->g2_values[i], &g2_affine);
    }

    /* It's the smallest power of 2 >= n1 */
    unsigned int max_scale = 0;
    while ((1ULL << max_scale) < n1)
        max_scale++;

    /* Initialize the KZGSettings struct */
    ret = c_kzg_malloc((void **)&out->fs, sizeof(FFTSettings));
    if (ret != C_KZG_OK) goto out_error;
    ret = new_fft_settings(out->fs, max_scale);
    if (ret != C_KZG_OK) goto out_error;
    ret = fft_g1(out->g1_values, g1_projective, true, n1, out->fs);
    if (ret != C_KZG_OK) goto out_error;
    ret = bit_reversal_permutation(out->g1_values, sizeof(g1_t), n1);
    if (ret != C_KZG_OK) goto out_error;

    goto out_success;

out_error:
    c_kzg_free(out->fs);
    c_kzg_free(out->g1_values);
    c_kzg_free(out->g2_values);
out_success:
    c_kzg_free(g1_projective);
    return ret;
}
*/
