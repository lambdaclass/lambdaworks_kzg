use crate::math::cyclic_group::IsGroup;
use crate::math::elliptic_curve::traits::FromAffine;
use crate::math::field::element::{FieldElement, LegendreSymbol};
use crate::math::field::extensions::quadratic::QuadraticExtensionFieldElement;
use crate::math::{errors::ByteConversionError, traits::ByteConversion};
use crate::G1Point;
use crate::MODULUS;
use crate::{BLS12381FieldElement, G2Point};
use std::cmp::Ordering;
use std::ops::Neg;

pub type QFE = FieldElement<crate::math::field::extensions::quadratic::QuadraticExtensionField<crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::field_extension::LevelOneResidue>>;

pub fn check_point_is_in_subgroup(point: &G1Point) -> bool {
    let inf = G1Point::neutral_element();
    let aux_point = point.operate_with_self(MODULUS);
    inf == aux_point
}

pub fn decompress_g1_point(input_bytes: &mut [u8; 48]) -> Result<G1Point, ByteConversionError> {
    let first_byte = input_bytes.first().unwrap();
    // We get the first 3 bits
    let prefix_bits = first_byte >> 5;
    let first_bit = (prefix_bits & 4_u8) >> 2;
    // If first bit is not 1, then the value is not compressed.
    if first_bit != 1 {
        return Err(ByteConversionError::ValueNotCompressed);
    }
    let second_bit = (prefix_bits & 2_u8) >> 1;
    // If the second bit is 1, then the compressed point is the
    // point at infinity and we return it directly.
    if second_bit == 1 {
        return Ok(G1Point::neutral_element());
    }
    let third_bit = prefix_bits & 1_u8;

    let first_byte_without_contorl_bits = (first_byte << 3) >> 3;
    input_bytes[0] = first_byte_without_contorl_bits;

    let x = BLS12381FieldElement::from_bytes_be(input_bytes)?;

    // We apply the elliptic curve formula to know the y^2 value.
    let y_squared = x.pow(3_u16) + BLS12381FieldElement::from(4);

    let (y_sqrt_1, y_sqrt_2) = &y_squared.sqrt().ok_or(ByteConversionError::InvalidValue)?;

    // we call "negative" to the greate root,
    // if the third bit is 1, we take this grater value.
    // Otherwise, we take the second one.
    let y = match (
        y_sqrt_1.representative().cmp(&y_sqrt_2.representative()),
        third_bit,
    ) {
        (Ordering::Greater, 0) => y_sqrt_2,
        (Ordering::Greater, _) => y_sqrt_1,
        (Ordering::Less, 0) => y_sqrt_1,
        (Ordering::Less, _) => y_sqrt_2,
        (Ordering::Equal, _) => y_sqrt_1,
    };

    let point =
        G1Point::from_affine(x, y.clone()).map_err(|_| ByteConversionError::InvalidValue)?;

    check_point_is_in_subgroup(&point)
        .then_some(point)
        .ok_or(ByteConversionError::PointNotInSubgroup)
}

pub fn decompress_g2_point(input_bytes: &mut [u8; 96]) -> Result<G2Point, ByteConversionError> {
    let binding = input_bytes[48..96].to_owned();
    let input0 = binding.as_slice();
    let input1 = &mut input_bytes[0..48];

    let first_byte = input1.first().unwrap();
    // We get the first 3 bits
    let prefix_bits = first_byte >> 5;
    let first_bit = (prefix_bits & 4_u8) >> 2;
    // If first bit is not 1, then the value is not compressed.
    if first_bit != 1 {
        return Err(ByteConversionError::ValueNotCompressed);
    }
    let second_bit = (prefix_bits & 2_u8) >> 1;
    // If the second bit is 1, then the compressed point is the
    // point at infinity and we return it directly.
    if second_bit == 1 {
        return Ok(G2Point::neutral_element());
    }
    let third_bit = prefix_bits & 1_u8;

    let first_byte_without_contorl_bits = (first_byte << 3) >> 3;
    input1[0] = first_byte_without_contorl_bits;

    let x0 = BLS12381FieldElement::from_bytes_be(input0).unwrap();
    let x1 = BLS12381FieldElement::from_bytes_be(input1).unwrap();
    let x: QFE = QuadraticExtensionFieldElement::new([x0, x1]);

    // TODO: calculate sqrt

    //let y0 = BLS12381FieldElement::from(2u64);
    //let y1 = BLS12381FieldElement::from(2u64);
    //let y = QuadraticExtensionFieldElement::new([y0, y1]);

    //let point = G2Point::from_affine(x, y).unwrap();
    todo!();
}

/// * `third_bit` - if 1, then the square root is the greater one, otherwise it is the smaller one.
pub fn sqrt_qfe(input: &QFE, third_bit: u8) -> Option<QFE> {
    // Algorithm 8, https://eprint.iacr.org/2012/685.pdf
    if *input == QFE::zero() {
        Some(QFE::zero())
    } else {
        let a = input.value()[0].clone();
        let b = input.value()[1].clone();
        if b == crate::BLS12381FieldElement::zero() {
            // second part is zero
            let (y_sqrt_1, y_sqrt_2) = a.sqrt()?;
            println!("y_sqrt_1: {:?}", y_sqrt_1);

            let y_aux = match (
                y_sqrt_1.representative().cmp(&y_sqrt_2.representative()),
                third_bit,
            ) {
                (Ordering::Greater, 0) => y_sqrt_2,
                (Ordering::Greater, _) => y_sqrt_1,
                (Ordering::Less, 0) => y_sqrt_1,
                (Ordering::Less, _) => y_sqrt_2,
                (Ordering::Equal, _) => y_sqrt_1,
            };

            Some(QFE::new([y_aux, crate::BLS12381FieldElement::zero()]))
        } else {
            // second part of the input field number is non-zero

            // instead of "sum" is -beta
            let alpha = a.pow(2u64) + b.pow(2u64);
            let gamma = alpha.legendre_symbol();
            match gamma {
                LegendreSymbol::One => {
                    let two = BLS12381FieldElement::from(2u64);
                    let two_inv = two.inv();
                    // calculate the square root of alpha
                    let (y_sqrt1, y_sqrt2) = alpha.sqrt()?;
                    println!("LegendreSymbol::One y_sqrt1: {:?}", y_sqrt1);

                    let mut delta = (a.clone() + y_sqrt1) * two_inv.clone();

                    let legendre_delta = delta.legendre_symbol();
                    if legendre_delta == LegendreSymbol::MinusOne {
                        delta = (a.clone() + y_sqrt2) * two_inv;
                    };
                    let (x_sqrt_1, x_sqrt_2) = delta.sqrt()?;
                    println!("LegendreSymbol::One x_sqrt_1: {:?}", x_sqrt_1);

                    let x_0 = match (
                        x_sqrt_1.representative().cmp(&x_sqrt_2.representative()),
                        third_bit,
                    ) {
                        (Ordering::Greater, 0) => x_sqrt_2,
                        (Ordering::Greater, _) => x_sqrt_1,
                        (Ordering::Less, 0) => x_sqrt_1,
                        (Ordering::Less, _) => x_sqrt_2,
                        (Ordering::Equal, _) => x_sqrt_1,
                    };
                    let x_1 = b * (two * x_0.clone()).inv();

                    Some(QFE::new([x_0, x_1]))
                }
                LegendreSymbol::MinusOne => {
                    println!("### LegendreSymbol::MinusOne");
                    None
                }
                LegendreSymbol::Zero => {
                    unreachable!("The input is zero, but we already handled this case.")
                }
            }
        }
    }
}

pub fn compress_g1_point(point: &G1Point) -> Result<[u8; 48], Vec<u8>> {
    let ret_vec = if *point == G1Point::neutral_element() {
        // point is at infinity
        let mut x_bytes = vec![0_u8; 48];
        x_bytes[0] |= 1 << 7;
        x_bytes[0] |= 1 << 6;
        x_bytes
    } else {
        // point is not at infinity
        let point_affine = point.to_affine();
        let x = point_affine.x();
        let y = point_affine.y();

        let mut x_bytes = x.to_bytes_be();

        // Set first bit to to 1 indicate this is compressed element.
        x_bytes[0] |= 1 << 7;

        let y_neg = y.neg();
        if y_neg.representative() < y.representative() {
            x_bytes[0] |= 1 << 5;
        }
        x_bytes
    };
    ret_vec.try_into()
}

#[cfg(test)]
mod tests {
    use crate::math::cyclic_group::IsGroup;
    use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::curve::BLS12381Curve;
    use crate::math::elliptic_curve::traits::{FromAffine, IsEllipticCurve};
    use crate::math::traits::ByteConversion;
    use crate::math::unsigned_integer::element::UnsignedInteger;
    use crate::{BLS12381FieldElement, G1Point};

    use super::{compress_g1_point, decompress_g1_point};

    #[test]
    fn test_zero_point() {
        let g1 = BLS12381Curve::generator();

        assert!(super::check_point_is_in_subgroup(&g1));
        let new_x = BLS12381FieldElement::zero();
        let new_y = BLS12381FieldElement::one() + BLS12381FieldElement::one();

        let false_point2 = G1Point::from_affine(new_x, new_y).unwrap();

        assert!(!super::check_point_is_in_subgroup(&false_point2));
    }

    #[test]
    fn test_g1_compress_generator() {
        let g = BLS12381Curve::generator();
        let mut compressed_g = compress_g1_point(&g).unwrap();
        let first_byte = compressed_g[0];

        let first_byte_without_control_bits = (first_byte << 3) >> 3;
        compressed_g[0] = first_byte_without_control_bits;

        let compressed_g_x = BLS12381FieldElement::from_bytes_be(&compressed_g).unwrap();
        let g_x = g.x();

        assert_eq!(*g_x, compressed_g_x);
    }

    #[test]
    fn test_g1_compress_point_at_inf() {
        let inf = G1Point::neutral_element();
        let compressed_inf = compress_g1_point(&inf).unwrap();
        let first_byte = compressed_inf[0];

        assert_eq!(first_byte >> 6, 3_u8);
    }

    #[test]
    fn test_compress_decompress_generator() {
        let g = BLS12381Curve::generator();
        let mut compressed_g = compress_g1_point(&g).unwrap();
        let decompressed_g = decompress_g1_point(&mut compressed_g).unwrap();

        assert_eq!(g, decompressed_g);
    }

    #[test]
    fn test_compress_decompress_2g() {
        let g = BLS12381Curve::generator();
        // calculate g point operate with itself
        let g_2 = g.operate_with_self(UnsignedInteger::<4>::from("2"));
        let mut compressed_g2 = compress_g1_point(&g_2).unwrap();
        let decompressed_g2 = decompress_g1_point(&mut compressed_g2).unwrap();

        assert_eq!(g_2, decompressed_g2);
    }

    #[test]
    fn test_sqrt_qfe() {
        let c1 = BLS12381FieldElement::from_hex(
            "0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e",
        ).unwrap();
        let c0 = BLS12381FieldElement::from_hex(
        "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8"
        ).unwrap();
        let qfe = super::QFE::new([c0, c1]);

        let b1 = BLS12381FieldElement::from_hex("0x4").unwrap();
        let b0 = BLS12381FieldElement::from_hex("0x4").unwrap();
        let qfe_b = super::QFE::new([b0, b1]);

        let cubic_value = qfe.pow(3_u64) + qfe_b;
        let root = super::sqrt_qfe(&cubic_value, 0).unwrap();

        let c0_expected = BLS12381FieldElement::from_hex("0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801").unwrap();
        let c1_expected = BLS12381FieldElement::from_hex("0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be").unwrap();
        let qfe_expected = super::QFE::new([c0_expected, c1_expected]);

        let value_root = root.value();
        let value_qfe_expected = qfe_expected.value();

        assert_eq!(value_root[0].clone(), value_qfe_expected[0].clone());
        assert_eq!(value_root[1].clone(), value_qfe_expected[1].clone());
    }

    #[test]
    fn test_sqrt_qfe_2() {
        let c0 = BLS12381FieldElement::from_hex("0x02").unwrap();
        let c1 = BLS12381FieldElement::from_hex("0x00").unwrap();
        let qfe = super::QFE::new([c0, c1]);

        let c0_expected = BLS12381FieldElement::from_hex("0x013a59858b6809fca4d9a3b6539246a70051a3c88899964a42bc9a69cf9acdd9dd387cfa9086b894185b9a46a402be73").unwrap();
        let c1_expected = BLS12381FieldElement::from_hex("0x02d27e0ec3356299a346a09ad7dc4ef68a483c3aed53f9139d2f929a3eecebf72082e5e58c6da24ee32e03040c406d4f").unwrap();
        let qfe_expected = super::QFE::new([c0_expected, c1_expected]);

        let b1 = BLS12381FieldElement::from_hex("0x4").unwrap();
        let b0 = BLS12381FieldElement::from_hex("0x4").unwrap();
        let qfe_b = super::QFE::new([b0, b1]);

        let root = super::sqrt_qfe(&(qfe.pow(3_u64) + qfe_b), 0).unwrap();

        let value_root = root.value();
        let value_qfe_expected = qfe_expected.value();

        assert_eq!(value_root[0].clone(), value_qfe_expected[0].clone());
        assert_eq!(value_root[1].clone(), value_qfe_expected[1].clone());
    }
}
