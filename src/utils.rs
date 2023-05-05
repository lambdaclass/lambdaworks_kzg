use crate::commitments::kzg::{FrElement, StructuredReferenceString, G1};
use crate::math::cyclic_group::IsGroup;
use crate::math::elliptic_curve::traits::IsEllipticCurve;
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
use crate::{BYTES_PER_BLOB, BYTES_PER_FIELD_ELEMENT};
use rand::Rng;

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
