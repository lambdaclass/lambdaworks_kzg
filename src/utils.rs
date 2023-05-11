use crate::commitments::kzg::{FrElement, StructuredReferenceString, G1};
use crate::compress::compress_g1_point;
use crate::math::cyclic_group::IsGroup;
use crate::math::elliptic_curve::traits::IsEllipticCurve;
use crate::math::errors::ByteConversionError;
use crate::math::unsigned_integer::element::U256;
use crate::math::{
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{
            curve::BLS12381Curve, pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
        },
        traits::IsPairing,
    },
    polynomial::Polynomial,
    traits::ByteConversion,
};
use crate::{
    G1Point, BYTES_PER_BLOB, BYTES_PER_FIELD_ELEMENT, FE, FIAT_SHAMIR_PROTOCOL_DOMAIN,
    FIELD_ELEMENTS_PER_BLOB,
};
use rand::Rng;

pub fn blob_to_polynomial(
    input_blob: &[u8; BYTES_PER_BLOB],
) -> Result<Polynomial<FrElement>, ByteConversionError>
where
    FrElement: ByteConversion,
{
    let mut coefficients = Vec::new();

    for elem_bytes in input_blob.chunks(BYTES_PER_FIELD_ELEMENT) {
        let f = FrElement::from_bytes_be(elem_bytes)?;
        coefficients.push(f);
    }

    Ok(Polynomial::new(&coefficients))
}

pub fn polynomial_to_blob(polynomial: &Polynomial<FrElement>) -> Vec<u8>
where
    FrElement: ByteConversion,
{
    let coefficients = polynomial.coefficients();

    coefficients
        .iter()
        .flat_map(|coef| coef.to_bytes_be())
        .collect()
}

pub fn polynomial_to_blob_with_size(
    polynomial: &Polynomial<FrElement>,
) -> Result<[u8; BYTES_PER_BLOB], Vec<u8>>
where
    FrElement: ByteConversion,
{
    let coefficients = polynomial.coefficients();
    let mut ret_vec: Vec<u8> = coefficients
        .iter()
        .flat_map(|coef| coef.to_bytes_be())
        .collect();

    let len = ret_vec.len();
    let remaining = (BYTES_PER_BLOB - len) / BYTES_PER_FIELD_ELEMENT;

    // pad with zeros until BYTES_PER_BLOB
    let zero = FE::zero();
    let zero_bytes = zero.to_bytes_be();

    for _i in 0..remaining {
        ret_vec.extend_from_slice(&zero_bytes);
    }
    ret_vec.try_into()
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

/// Return the Fiat-Shamir challenge required to verify `blob` and
/// `commitment`.
///
/// # Params
///
/// - `blob` - A blob
/// - `commitment` - A commitment
///
/// # Returns
///
/// FrElement corresponding to the field element value.
pub fn compute_challenge(
    blob: &[u8; BYTES_PER_BLOB],
    commitment_g1: &G1Point,
) -> Result<FrElement, Vec<u8>> {
    // insert the values in the string, hash and get the field element
    // concat:
    // - FIAT_SHAMIR_PROTOCOL_DOMAIN
    // - FIELD_ELEMENTS_PER_BLOB as litlle-endian number
    // - 0 as u64
    // - blob (this is BYTES_PER_BLOB bytes)
    // - g1 point

    let input_hash = FIAT_SHAMIR_PROTOCOL_DOMAIN
        .into_iter()
        .chain(FIELD_ELEMENTS_PER_BLOB.to_le_bytes().into_iter())
        .chain(0_u64.to_le_bytes().into_iter())
        .chain(blob.iter().copied())
        .chain(compress_g1_point(commitment_g1)?.into_iter())
        .collect::<Vec<u8>>();
    hash_field_unsafe(&input_hash)
}

/// Hashes the input sting and returns the field element corresponding to
/// the hash coverted to field
fn hash_field_unsafe(input_slice: &[u8]) -> Result<FrElement, Vec<u8>> {
    let ret_hash = sha256::digest(input_slice);
    let mut bytes_hash = [0u8; 32];
    hex::decode_to_slice(&ret_hash, &mut bytes_hash as &mut [u8]).map_err(|_| Vec::new())?;
    // FIXME! This should be changed to a hash to field method
    FrElement::from_bytes_be(&bytes_hash).map_err(|_| Vec::new())
}

#[cfg(test)]
mod tests {
    use super::{blob_to_polynomial, polynomial_to_blob_with_size};
    use crate::commitments::kzg::FrElement;
    use crate::math::field::element::FieldElement;
    use crate::math::polynomial::Polynomial;
    use crate::FE;

    #[test]
    fn test_poly_to_blob_and_viceversa() {
        let polynomial = Polynomial::<FrElement>::new(&[FieldElement::one()]);
        let blob = polynomial_to_blob_with_size(&polynomial).unwrap();
        let poly_from_blob = blob_to_polynomial(&blob).unwrap();

        let one = FE::from(1);
        let y_1 = polynomial.evaluate(&one);
        let y_2 = poly_from_blob.evaluate(&one);

        assert_eq!(y_1, y_2);
    }
}
