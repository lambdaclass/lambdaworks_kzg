use crate::compress::{decompress_g1_point, decompress_g2_point};
use crate::{blst_fp, blst_fp2, blst_p1, blst_p2, G2Point, KZGSettings, G1};
use core::ptr::null_mut;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines};
use std::marker;
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
pub fn load_trusted_setup_file(path: &str) -> io::Result<KZGSettings> {
    let mut lines = read_lines(path)?;

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
        ))??
        .parse::<usize>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    let num_g2_points = lines
        .next()
        .ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format",
        ))??
        .parse::<usize>()
        .map_err(|_| std::io::ErrorKind::InvalidData)?;

    println!("num_g1_points: {num_g1_points}, num_g2_points: {num_g2_points}");

    let num_total_points = num_g1_points + num_g2_points;

    // read all g1 points
    for (pos, line) in lines.enumerate() {
        if pos < num_g1_points {
            // read g1 point
            let Ok(line_string) = line  else {
                return Err(std::io::ErrorKind::InvalidData.into());
            };
            hex::decode_to_slice(line_string, &mut g1_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            let g1_point =
                decompress_g1_point(&mut g1_bytes).map_err(|_| std::io::ErrorKind::InvalidData)?;

            g1_points.push(g1_point);
        } else if pos < num_total_points {
            // read g2 point
            let Ok(line_string) = line  else {
                return Err(std::io::ErrorKind::InvalidData.into());
            };
            hex::decode_to_slice(line_string, &mut g2_bytes)
                .map_err(|_| std::io::ErrorKind::InvalidData)?;

            let g2_point =
                decompress_g2_point(&mut g2_bytes).map_err(|_| std::io::ErrorKind::InvalidData)?;

            g2_points.push(g2_point);
        } else {
            // all the points were already parsed
            break;
        }
    }

    let mut g1_values_vec: Vec<blst_p1> = g1_points
        .iter()
        .map(|v| blst_p1 {
            x: blst_fp {
                l: v.x().representative().limbs,
            },
            y: blst_fp {
                l: v.y().representative().limbs,
            },
            z: blst_fp {
                l: v.z().representative().limbs,
            },
        })
        .collect();
    let g1_values = g1_values_vec.as_mut_ptr();
    std::mem::forget(g1_values_vec);

    let mut g2_values_vec: Vec<blst_p2> = g2_points
        .iter()
        .map(|v| {
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
        })
        .collect();

    let g2_values = g2_values_vec.as_mut_ptr();
    std::mem::forget(g2_values_vec);

    let settings = KZGSettings {
        fs: null_mut(),
        g1_values,
        g2_values,
        _marker: marker::PhantomData,
        _marker2: marker::PhantomData,
        _marker3: marker::PhantomData,
    };

    let len_vec_g1 = g1_points.len();
    let len_vec_g2 = g2_points.len();
    println!("len_vec_g1: {len_vec_g1}, len_vec_g2: {len_vec_g2}");

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

#[cfg(test)]
mod tests {
    use super::load_trusted_setup_file;

    #[test]
    fn test_read_srs() {
        load_trusted_setup_file("test/trusted_setup_4.txt").unwrap();
    }
}
