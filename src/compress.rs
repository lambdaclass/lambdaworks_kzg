use crate::math::cyclic_group::IsGroup;
use crate::math::elliptic_curve::traits::FromAffine;
use crate::math::field::extensions::quadratic::QuadraticExtensionFieldElement;
use crate::math::{errors::ByteConversionError, traits::ByteConversion};
use crate::sqrt::select_sqrt_value_from_third_bit;
use crate::BLS12381TwistCurveFieldElement;
use crate::G1Point;
use crate::MODULUS;
use crate::{BLS12381FieldElement, G2Point};
use std::ops::Neg;

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
    let y = select_sqrt_value_from_third_bit(y_sqrt_1.clone(), y_sqrt_2.clone(), third_bit);
    let point = G1Point::from_affine(x, y).map_err(|_| ByteConversionError::InvalidValue)?;

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

    let first_byte_without_contorl_bits = (first_byte << 3) >> 3;
    input1[0] = first_byte_without_contorl_bits;

    let x0 = BLS12381FieldElement::from_bytes_be(input0).unwrap();
    let x1 = BLS12381FieldElement::from_bytes_be(input1).unwrap();
    let x: BLS12381TwistCurveFieldElement = QuadraticExtensionFieldElement::new([x0, x1]);

    let b_param_qfe = get_four_qfe();
    let y = super::sqrt::sqrt_qfe(&(x.pow(3_u64) + b_param_qfe), 0)
        .ok_or(ByteConversionError::InvalidValue)?;
    G2Point::from_affine(x, y).map_err(|_| ByteConversionError::InvalidValue)
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

fn get_four_qfe() -> BLS12381TwistCurveFieldElement {
    let b1 = BLS12381FieldElement::from_hex("0x4").unwrap();
    let b0 = BLS12381FieldElement::from_hex("0x4").unwrap();
    super::BLS12381TwistCurveFieldElement::new([b0, b1])
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
    fn short_test_compress_and_decompress_point() {
        let line = "8d0c6eeadd3f8529d67246f77404a4ac2d9d7fd7d50cf103d3e6abb9003e5e36d8f322663ebced6707a7f46d97b7566d";
        let bytes = hex::decode(line).unwrap();
        let mut input_bytes: [u8; 48] = bytes.try_into().unwrap();
        let point = decompress_g1_point(&mut input_bytes).unwrap();
        let compressed = compress_g1_point(&point).unwrap();
        let hex_string = hex::encode(compressed);

        assert_eq!("8d0c6eeadd3f8529d67246f77404a4ac2d9d7fd7d50cf103d3e6abb9003e5e36d8f322663ebced6707a7f46d97b7566d", &hex_string);
    }
}
