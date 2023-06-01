use crate::FieldElement;
use crate::{
    G1Point, G2Point, BYTES_PER_BLOB, BYTES_PER_FIELD_ELEMENT, FE, FIAT_SHAMIR_PROTOCOL_DOMAIN,
    FIELD_ELEMENTS_PER_BLOB, G1, RANDOM_CHALLENGE_KZG_BATCH_DOMAIN,
};
use itertools::izip;
use lambdaworks_crypto::commitments::kzg::StructuredReferenceString;
use lambdaworks_math::cyclic_group::IsGroup;
use lambdaworks_math::elliptic_curve::short_weierstrass::curves::bls12_381::field_extension::LevelThreeResidue;
use lambdaworks_math::elliptic_curve::traits::{Compress, IsEllipticCurve};
use lambdaworks_math::errors::ByteConversionError;
use lambdaworks_math::field::extensions::quadratic::QuadraticExtensionField;
use lambdaworks_math::unsigned_integer::element::U256;
use lambdaworks_math::{
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{
            curve::BLS12381Curve, pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
        },
        traits::IsPairing,
    },
    polynomial::Polynomial,
    traits::ByteConversion,
};
use rand::Rng;

pub fn blob_to_polynomial(
    input_blob: &[u8; BYTES_PER_BLOB],
) -> Result<Polynomial<FE>, ByteConversionError>
where
    FE: ByteConversion,
{
    let mut coefficients = Vec::new();

    for elem_bytes in input_blob.chunks(BYTES_PER_FIELD_ELEMENT) {
        let f = FE::from_bytes_be(elem_bytes)?;
        coefficients.push(f);
    }

    Ok(Polynomial::new(&coefficients))
}

#[must_use]
pub fn polynomial_to_blob(polynomial: &Polynomial<FE>) -> Vec<u8>
where
    FE: ByteConversion,
{
    let coefficients = polynomial.coefficients();

    coefficients
        .iter()
        .flat_map(ByteConversion::to_bytes_be)
        .collect()
}

pub fn polynomial_to_blob_with_size(
    polynomial: &Polynomial<FE>,
) -> Result<[u8; BYTES_PER_BLOB], Vec<u8>>
where
    FE: ByteConversion,
{
    let coefficients = polynomial.coefficients();
    let mut ret_vec: Vec<u8> = coefficients
        .iter()
        .flat_map(ByteConversion::to_bytes_be)
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
#[must_use]
pub fn create_srs() -> StructuredReferenceString<
    <BLS12381AtePairing as IsPairing>::G1Point,
    <BLS12381AtePairing as IsPairing>::G2Point,
> {
    let mut rng = rand::thread_rng();
    let toxic_waste = FE::new(U256 {
        limbs: [
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
        ],
    });
    let g1 = BLS12381Curve::generator();
    let g2 = BLS12381TwistCurve::generator();
    let powers_main_group: Vec<G1> = (0..100_u128)
        .map(|exponent| g1.operate_with_self(toxic_waste.pow(exponent).representative()))
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
/// FE corresponding to the field element value.
pub fn compute_challenge(
    blob: &[u8; BYTES_PER_BLOB],
    commitment_g1: &G1Point,
) -> Result<FE, Vec<u8>> {
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
        .chain(
            G1Point::compress_g1_point(commitment_g1)
                .map_err(|_| Vec::new())?
                .into_iter(),
        )
        .collect::<Vec<u8>>();
    hash_field_unsafe(&input_hash)
}

/// Hashes the input sting and returns the field element corresponding to
/// the hash coverted to field
fn hash_field_unsafe(input_slice: &[u8]) -> Result<FE, Vec<u8>> {
    let ret_hash = sha256::digest(input_slice);
    let mut bytes_hash = [0u8; 32];
    hex::decode_to_slice(&ret_hash, &mut bytes_hash as &mut [u8]).map_err(|_| Vec::new())?;
    // FIXME! This should be changed to a hash to field method
    FE::from_bytes_be(&bytes_hash).map_err(|_| Vec::new())
}

fn compute_powers(x: &FE, n: usize) -> Vec<FE> {
    let mut current_power = FE::one();
    let mut powers = Vec::with_capacity(n);
    for _i in 0..n {
        powers.push(current_power.clone());
        current_power = current_power * x;
    }
    powers
}

pub fn compute_r_powers(
    commitments_g1: &[G1Point],
    zs_fr: &[FE],
    ys_fr: &[FE],
    proofs_g1: &[G1Point],
) -> Result<Vec<FE>, Vec<u8>> {
    let n = commitments_g1.len();
    let mut bytes: Vec<u8> = RANDOM_CHALLENGE_KZG_BATCH_DOMAIN
        .into_iter()
        .chain(FIELD_ELEMENTS_PER_BLOB.to_le_bytes().into_iter())
        .chain(n.to_le_bytes().into_iter())
        .collect();

    for (commitment, z, y, proof) in izip!(
        commitments_g1.iter(),
        zs_fr.iter(),
        ys_fr.iter(),
        proofs_g1.iter()
    ) {
        // TODO: make it beautiful
        let mut input_hash = Vec::<u8>::new();
        input_hash.extend_from_slice(
            G1Point::compress_g1_point(commitment)
                .map_err(|_| Vec::new())?
                .as_slice(),
        );
        input_hash.extend_from_slice(&z.to_bytes_be());
        input_hash.extend_from_slice(&y.to_bytes_be());
        input_hash.extend_from_slice(
            G1Point::compress_g1_point(proof)
                .map_err(|_| Vec::new())?
                .as_slice(),
        );

        bytes.append(&mut input_hash);
    }

    // Now let's create the challenge!
    let hash_fr = hash_field_unsafe(&bytes)?;
    Ok(compute_powers(&hash_fr, n))
}

/// Perform pairings and test whether the outcomes are equal in `G_T`.
///
/// Tests whether `e(a1, a2) == e(b1, b2)`.
///
/// # Params
///
/// * `a1` - A G1 group point for the first pairing
/// * `a2` - A G2 group point for the first pairing
/// * `b1` - A G1 group point for the second pairing
/// * `b2` - A G2 group point for the second pairing
///
/// # Returns
///
/// * true  The pairings were equal
/// * false The pairings were not equal
#[must_use]
pub fn pairings_verify(a1: &G1Point, a2: &G2Point, b1: &G1Point, b2: &G2Point) -> bool {
    // As an optimisation, we want to invert one of the pairings,
    // so we negate one of the points.
    let a1neg = a1.neg();
    let aa1 = a1neg.to_affine();
    let bb1 = b1.to_affine();
    let aa2 = a2.to_affine();
    let bb2 = b2.to_affine();

    let pairs = vec![(&aa1, &aa2), (&bb1, &bb2)];
    let gt_point = BLS12381AtePairing::compute_batch(&pairs);
    FieldElement::<QuadraticExtensionField<LevelThreeResidue>>::one() == gt_point
}

#[cfg(test)]
mod tests {
    use super::{blob_to_polynomial, polynomial_to_blob_with_size};
    use crate::FE;
    use lambdaworks_math::field::element::FieldElement;
    use lambdaworks_math::polynomial::Polynomial;

    #[test]
    fn test_poly_to_blob_and_viceversa() {
        let polynomial = Polynomial::<FE>::new(&[FieldElement::one()]);
        let blob = polynomial_to_blob_with_size(&polynomial).unwrap();
        let poly_from_blob = blob_to_polynomial(&blob).unwrap();

        let one = FE::from(1);
        let y_1 = polynomial.evaluate(&one);
        let y_2 = poly_from_blob.evaluate(&one);

        assert_eq!(y_1, y_2);
    }
}
