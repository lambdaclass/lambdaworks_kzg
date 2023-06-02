use crate::{
    blst_fp, blst_fp2, blst_p1, blst_p2, traits::Compress, BLS12381FieldElement,
    BLS12381TwistCurveFieldElement, G1Point, G2Point, KZGSettings, G1, NUM_G1_POINTS,
    NUM_G2_POINTS,
};
use core::ptr::null_mut;
use lambdaworks_crypto::commitments::kzg::StructuredReferenceString;
use lambdaworks_math::cyclic_group::IsGroup;
use lambdaworks_math::errors::ByteConversionError;
use lambdaworks_math::{elliptic_curve::traits::FromAffine, traits::ByteConversion};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Helper function that reads a file line by line and
/// returns an iterator over the lines.
pub fn read_lines<P>(path: P) -> io::Result<std::io::Lines<BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(path)?;
    Ok(BufReader::new(file).lines())
}

pub fn load_trusted_setup_file_to_g1_points_and_g2_points(
    mut lines: std::str::Lines,
) -> io::Result<(Vec<G1>, Vec<G2Point>)> {
    let mut g1_bytes: [u8; crate::BYTES_PER_G1_POINT] = [0; crate::BYTES_PER_G1_POINT];
    let mut g1_points: Vec<G1> = Vec::new();

    let mut g2_bytes: [u8; crate::BYTES_PER_G2_POINT] = [0; crate::BYTES_PER_G2_POINT];
    let mut g2_points: Vec<G2Point> = Vec::new();

    // Read the number of g1 points
    let num_g1_points = lines
        .next()
        .ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format",
        ))?
        .parse::<usize>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    let num_g2_points = lines
        .next()
        .ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format",
        ))?
        .parse::<usize>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    let num_total_points = num_g1_points + num_g2_points;

    // read all g1 points
    for (pos, line) in lines.enumerate() {
        if pos < num_g1_points {
            // read g1 point
            hex::decode_to_slice(line, &mut g1_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            let g1_point = G1Point::decompress_g1_point(&mut g1_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            g1_points.push(g1_point);
        } else if pos < num_total_points {
            // read g2 point
            hex::decode_to_slice(line, &mut g2_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            let g2_point = G1Point::decompress_g2_point(&mut g2_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            g2_points.push(g2_point);
        } else {
            // all the points were already parsed
            break;
        }
    }

    Ok((g1_points, g2_points))
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
pub fn load_trusted_setup_file(lines: std::str::Lines<'_>) -> io::Result<KZGSettings> {
    let (g1_points, g2_points) = load_trusted_setup_file_to_g1_points_and_g2_points(lines)?;

    let mut g1_values_vec: Vec<blst_p1> = g1_points.iter().map(g1_point_to_blst_p1).collect();
    let g1_values = g1_values_vec.as_mut_ptr();
    std::mem::forget(g1_values_vec);

    let mut g2_values_vec: Vec<blst_p2> = g2_points.iter().map(g2_point_to_blst_p2).collect();

    let g2_values = g2_values_vec.as_mut_ptr();
    std::mem::forget(g2_values_vec);

    let settings = KZGSettings {
        fs: null_mut(),
        g1_values,
        g2_values,
    };

    /* TODO: add this to the KZGSettings struct
    ret = new_fft_settings(out->fs, max_scale);
    if (ret != C_KZG_OK) goto out_error;
    ret = fft_g1(out->g1_values, g1_projective, true, n1, out->fs);
    if (ret != C_KZG_OK) goto out_error;
    ret = bit_reversal_permutation(out->g1_values, sizeof(g1_t), n1);
    if (ret != C_KZG_OK) goto out_error;
    */
    // TODO - return the KZGSettings:
    // convert vector of g1 and g2 point to lambdaworks format
    Ok(settings)
}

#[must_use]
pub fn g1_point_to_blst_p1(v: &G1) -> blst_p1 {
    if v.is_neutral_element() {
        return blst_p1 {
            x: blst_fp { l: [0; 6] },
            y: blst_fp { l: [0; 6] },
            z: blst_fp {
                l: [0, 0, 0, 0, 0, 1],
            },
        };
    }

    blst_p1 {
        x: blst_fp {
            l: v.x().representative().limbs,
        },
        y: blst_fp {
            l: v.y().representative().limbs,
        },
        z: blst_fp {
            l: v.z().representative().limbs,
        },
    }
}

pub fn blst_p1_to_g1_point(point: &blst_p1) -> Result<G1, ByteConversionError> {
    let x_field = point
        .x
        .l
        .iter()
        .flat_map(|e| e.to_be_bytes())
        .collect::<Vec<u8>>();
    let x = BLS12381FieldElement::from_bytes_be(&x_field)?;

    let y_field = point
        .y
        .l
        .iter()
        .flat_map(|e| e.to_be_bytes())
        .collect::<Vec<u8>>();
    let y = BLS12381FieldElement::from_bytes_be(&y_field)?;
    G1::from_affine(x, y).map_err(|_| ByteConversionError::FromBEBytesError)
}

#[must_use]
pub fn g2_point_to_blst_p2(v: &G2Point) -> blst_p2 {
    let vx = v.to_affine().x().value().clone();
    let x = blst_fp2 {
        fp: [
            blst_fp {
                l: vx[0].representative().limbs,
            },
            blst_fp {
                l: vx[1].representative().limbs,
            },
        ],
    };

    let vy = v.to_affine().y().value().clone();
    let y = blst_fp2 {
        fp: [
            blst_fp {
                l: vy[0].representative().limbs,
            },
            blst_fp {
                l: vy[1].representative().limbs,
            },
        ],
    };

    let vz = v.to_affine().z().value().clone();
    let z = blst_fp2 {
        fp: [
            blst_fp {
                l: vz[0].representative().limbs,
            },
            blst_fp {
                l: vz[1].representative().limbs,
            },
        ],
    };

    blst_p2 { x, y, z }
}

pub fn blst_p2_to_g2_point(point: &blst_p2) -> Result<G2Point, ByteConversionError> {
    let [x0, x1] = point.x.fp;
    let [y0, y1] = point.y.fp;

    let x0_field = BLS12381FieldElement::from_bytes_be(
        &x0.l
            .iter()
            .flat_map(|e| e.to_be_bytes())
            .collect::<Vec<u8>>(),
    )?;
    let x1_field = BLS12381FieldElement::from_bytes_be(
        &x1.l
            .iter()
            .flat_map(|e| e.to_be_bytes())
            .collect::<Vec<u8>>(),
    )?;
    let y0_field = BLS12381FieldElement::from_bytes_be(
        &y0.l
            .iter()
            .flat_map(|e| e.to_be_bytes())
            .collect::<Vec<u8>>(),
    )?;
    let y1_field = BLS12381FieldElement::from_bytes_be(
        &y1.l
            .iter()
            .flat_map(|e| e.to_be_bytes())
            .collect::<Vec<u8>>(),
    )?;

    let x = BLS12381TwistCurveFieldElement::new([x0_field, x1_field]);
    let y = BLS12381TwistCurveFieldElement::new([y0_field, y1_field]);
    G2Point::from_affine(x, y).map_err(|_| ByteConversionError::FromBEBytesError)
}

#[must_use]
pub fn vecs_to_structured_reference_string(
    g1_points: &[G1],
    g2_points: &[G2Point],
) -> StructuredReferenceString<G1, G2Point> {
    let g2_points_arr = [g2_points[0].clone(), g2_points[1].clone()];
    StructuredReferenceString::<G1, G2Point>::new(g1_points, &g2_points_arr)
}

pub fn kzgsettings_to_structured_reference_string(
    s: &KZGSettings,
) -> Result<StructuredReferenceString<G1, G2Point>, ByteConversionError> {
    let g1_points_slice: [blst_p1; NUM_G1_POINTS] = unsafe { *s.g1_values.cast() };
    let g2_points_slice: [blst_p2; NUM_G2_POINTS] = unsafe { *s.g2_values.cast() };

    // construct g1_points vec
    let g1_points: Result<Vec<G1>, ByteConversionError> =
        g1_points_slice.iter().map(blst_p1_to_g1_point).collect();
    let g1_points = g1_points?;

    // construct gs_points vec
    let g2_points: Result<Vec<G2Point>, ByteConversionError> =
        g2_points_slice.iter().map(blst_p2_to_g2_point).collect();
    let g2_points = g2_points?;

    let g2_points_arr = [g2_points[0].clone(), g2_points[1].clone()];

    Ok(StructuredReferenceString::<G1, G2Point>::new(
        &g1_points,
        &g2_points_arr,
    ))
}

#[cfg(test)]
mod tests {
    use super::load_trusted_setup_file;

    #[test]
    fn test_read_srs() {
        let lines = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/trusted_setup_4.txt"
        ));
        let lines = lines.lines();
        load_trusted_setup_file(lines).unwrap();
    }
}
