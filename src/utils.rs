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

pub fn get_point_from_bytes(input_bytes: &mut [u8; 48]) -> Result<G1Point, ByteConversionError> {
    let first_byte = input_bytes.first().unwrap();

    // We get the first 3 bits t
    let prefix_bits = first_byte >> 5;

    // let _first_bit = prefix_bits & 4_u8;
    let second_bit = prefix_bits & 2_u8;
    let third_bit = prefix_bits & 1_u8;

    if second_bit == 1 {
        return Ok(G1Point::neutral_element());
    }
    let first_byte_without_contorl_bits = (first_byte << 3) >> 3;
    input_bytes[0] = first_byte_without_contorl_bits;

    let x = BLS12381FieldElement::from_bytes_be(input_bytes)?;

    // We apply the elliptic curve formula to know the y^2 value.
    let y_squared = x.pow(3_u16) + BLS12381FieldElement::from(4);
    let y = sqrt_fr(&y_squared);

    let point = G1Point::from_affine(x, y).map_err(|_| ByteConversionError::InvalidValue)?;

    if check_point_is_in_subgroup(&point) {
        Ok(point)
    } else {
        Err(ByteConversionError::PointNotInSubgroup)
    }

    // * sacar los 3 bits mas significativos y guardarlos
    // el segundo indica si es el punto en el infinito
    // el terceer bit más significativo sirve para la raiz cuadrada
    // armo el field element from bytes big endian -> x

    // con x, hago la raíz cuadrada de (x^3 + 4)
    // tengo 2 raíces cuadradas, el 3er bit mas significativo,
    // me dice cuál obtener

    // creo el punto nuevo con el método g1_point = new_affine(x, y);

    // DONE!!!
    // chequear que el punto esté en el subgrupo:
    // se lo hace multiplicando el punto por el valor del módulo
    // si da el punto en el infinito, pertenece al subgrupo
    //

    // sacar los 3 bits
    // let a = BLS12381FieldElement::from_bytes_be(&input_bytes);
}

pub fn sqrt_fr(field: &BLS12381FieldElement) -> BLS12381FieldElement {
    field.inv().pow(2_u32)
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
    use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::curve::BLS12381Curve;
    use crate::math::elliptic_curve::traits::{FromAffine, IsEllipticCurve};
    use crate::{BLS12381FieldElement, G1Point};

    #[test]
    fn test_zero_point() {
        let g1 = BLS12381Curve::generator();

        assert!(super::check_point_is_in_subgroup(&g1));
        let new_x = BLS12381FieldElement::zero();
        let new_y = BLS12381FieldElement::one() + BLS12381FieldElement::one();

        // y2=x3+4

        let false_point2 = G1Point::from_affine(new_x, new_y).unwrap();

        assert!(!super::check_point_is_in_subgroup(&false_point2));
    }
}
