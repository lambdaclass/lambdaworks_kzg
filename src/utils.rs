use std::ops::Neg;

use crate::commitments::kzg::{FrElement, StructuredReferenceString, G1};
use crate::math::cyclic_group::IsGroup;
use crate::math::elliptic_curve::traits::{FromAffine, IsEllipticCurve};
use crate::math::unsigned_integer::element::U256;
use crate::math::{
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{
            curve::BLS12381Curve, pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
        },
        traits::IsPairing,
    },
    errors::ByteConversionError,
    polynomial::Polynomial,
    traits::ByteConversion,
};
use crate::G1Point;
use crate::{BLS12381FieldElement, BYTES_PER_BLOB, BYTES_PER_FIELD_ELEMENT};
use rand::Rng;

const MODULUS: U256 =
    U256::from("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");

pub fn check_point_is_in_subgroup(point: &G1Point) -> bool {
    let inf = G1Point::neutral_element();
    let aux_point = point.operate_with_self(MODULUS);
    inf == aux_point
}

#[allow(dead_code)]
pub fn blob_to_polynomial(
    input_blob: &[u8; BYTES_PER_BLOB],
) -> Result<Polynomial<FrElement>, crate::math::errors::ByteConversionError>
where
    FrElement: ByteConversion,
{
    let mut coefficients = Vec::new();

    for elem_bytes in input_blob.chunks(BYTES_PER_FIELD_ELEMENT) {
        let f = FrElement::from_bytes_le(elem_bytes)?;
        coefficients.push(f);
    }

    Ok(Polynomial::new(&coefficients))
}

pub fn decompress_g1_point(input_bytes: &mut [u8; 48]) -> Result<G1Point, ByteConversionError> {
    let first_byte = input_bytes.first().unwrap();

    // We get the first 3 bits
    let prefix_bits = first_byte >> 5;

    let first_bit = prefix_bits & 4_u8;
    // If first bit is not 1, then the value is not compressed.
    if first_bit != 1 {
        return Err(ByteConversionError::InvalidValue);
    }
    let second_bit = prefix_bits & 2_u8;
    let third_bit = prefix_bits & 1_u8;

    // If the second bit is 1, then the compressed point is the
    // point at infinity and we return it directly.
    if second_bit == 1 {
        return Ok(G1Point::neutral_element());
    }
    let first_byte_without_contorl_bits = (first_byte << 3) >> 3;
    input_bytes[0] = first_byte_without_contorl_bits;

    let x = BLS12381FieldElement::from_bytes_be(input_bytes)?;

    // We apply the elliptic curve formula to know the y^2 value.
    let y_squared = x.pow(3_u16) + BLS12381FieldElement::from(4);
    let y_abs = sqrt_fr(&y_squared);

    let y = if third_bit == 0 { y_abs } else { -y_abs };

    let point = G1Point::from_affine(x, y).map_err(|_| ByteConversionError::InvalidValue)?;

    if check_point_is_in_subgroup(&point) {
        Ok(point)
    } else {
        Err(ByteConversionError::PointNotInSubgroup)
    }
}

pub fn sqrt_fr(field: &BLS12381FieldElement) -> BLS12381FieldElement {
    field.inv().pow(2_u32)
}

pub fn compress_g1_point(point: &G1Point) -> Result<Vec<u8>, ByteConversionError> {
    let x = point.x();
    let y = point.y();

    let mut x_bytes = x.to_bytes_be();

    // Set first bit to to 1 indicate this is compressed element.
    x_bytes[0] |= 1 << 7;

    // Check that it is the point at infinity, if it is return
    if *point == G1Point::neutral_element() {
        x_bytes[0] |= 1 << 6;
        return Ok(x_bytes);
    };

    let y_neg = y.neg();
    if y_neg.representative() < y.representative() {
        x_bytes[0] |= 1 << 5;
    }
    Ok(x_bytes)
}

/// Helper function to create SRS. Once the deserialization of
/// the SRS is done, this function should be removed.
pub fn create_srs() -> StructuredReferenceString<
    <BLS12381AtePairing as IsPairing>::G1Point,
    <BLS12381AtePairing as IsPairing>::G2Point,
> {
    let mut rng = rand::thread_rng();
    let toxic_waste = FrElement::new(U256 {
        limbs: [
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
        ],
    });
    let g1 = BLS12381Curve::generator();
    let g2 = BLS12381TwistCurve::generator();
    let powers_main_group: Vec<G1> = (0..100)
        .map(|exponent| g1.operate_with_self(toxic_waste.pow(exponent as u128).representative()))
        .collect();
    let powers_secondary_group = [
        g2.clone(),
        g2.operate_with_self(toxic_waste.representative()),
    ];
    StructuredReferenceString::new(&powers_main_group, &powers_secondary_group)
}

#[cfg(test)]
mod tests {
    use crate::math::cyclic_group::IsGroup;
    use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::curve::BLS12381Curve;
    use crate::math::elliptic_curve::traits::{FromAffine, IsEllipticCurve};
    use crate::math::traits::ByteConversion;
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
        let first_byte = compressed_g.first().unwrap();

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
        let first_byte = compressed_inf.first().unwrap();

        assert_eq!(*first_byte >> 6, 3_u8);
    }

    #[test]
    fn test_compress_decompress_generator() {
        let g = BLS12381Curve::generator();
        let compressed_g = compress_g1_point(&g).unwrap();
        let mut compressed_g_slice: [u8; 48] = compressed_g.try_into().unwrap();

        let decompressed_g = decompress_g1_point(&mut compressed_g_slice).unwrap();

        assert_eq!(g, decompressed_g);
    }
}
