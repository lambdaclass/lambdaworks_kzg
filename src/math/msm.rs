use std::mem::size_of;

use crate::{
    commitments::kzg::{FrElement, FrField, G1},
    math::{cyclic_group::IsGroup, field::traits::IsField},
};

const WINDOW_SIZE: usize = 3;

/// This function computes the multiscalar multiplication (MSM).
///
/// Assume a group G of order r is given.
/// Let `hidings = [g_1, ..., g_n]` be a tuple of group points in G and
/// let `cs = [k_1, ..., k_n]` be a tuple of scalars in the Galois field GF(r).
///
/// Then, with additive notation, `msm(cs, hidings)` computes k_1 * g_1 + .... + k_n * g_n.
///
/// If `hidings` and `cs` are empty, then `msm` returns the zero element of the group.
///
/// Panics if `cs` and `hidings` have different lengths.
pub fn msm(cs: &[FrElement], hidings: &[G1]) -> G1 {
    debug_assert_eq!(
        cs.len(),
        hidings.len(),
        "Slices `cs` and `hidings` must be of the same length to compute `msm`."
    );

    let num_windows_per_limb = size_of::<<FrField as IsField>::BaseType>() / WINDOW_SIZE;
    let mut buckets = vec![G1::neutral_element(); (1 << WINDOW_SIZE) - 1];

    cs.iter()
        .zip(hidings)
        .map(|(c, hiding)| {
            c.representative().limbs.iter().for_each(|limb| {
                (0..num_windows_per_limb).rev().for_each(|window_idx| {
                    // Put in the right bucket the corresponding ps[i] for the current window.

                    let m_ij =
                        ((limb >> (window_idx * WINDOW_SIZE)) & ((1 << WINDOW_SIZE) - 1)) as usize;
                    if m_ij != 0 {
                        buckets[m_ij - 1] = buckets[m_ij - 1].operate_with(hiding);
                    }
                });
            });

            // Do the reduction step for the buckets.
            buckets
                .iter_mut()
                .rev()
                .scan(G1::neutral_element(), |m, b| {
                    *m = m.operate_with(b); // Reduction step.
                    *b = G1::neutral_element(); // Cleanup bucket slot to reuse in the next window.
                    Some(m.clone())
                })
                .reduce(|g, m| g.operate_with(&m))
                .unwrap_or_else(G1::neutral_element)
        })
        .reduce(|t, g| t.operate_with_self(1_u64 << WINDOW_SIZE).operate_with(&g))
        .unwrap_or_else(G1::neutral_element)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::curve::BLS12381Curve;
    use crate::math::elliptic_curve::short_weierstrass::point::ShortWeierstrassProjectivePoint;
    use crate::math::elliptic_curve::traits::IsEllipticCurve;
    use crate::math::field::fields::u64_prime_field::U64FieldElement;

    const ORDER_R: u64 = 5;
    type FE = U64FieldElement<ORDER_R>;

    #[test]
    fn msm_11_is_1_over_elliptic_curves() {
        let c: [FrElement; 1] = [FrElement::from(1)];
        let hiding = [BLS12381Curve::generator()];
        assert_eq!(msm(&c, &hiding), BLS12381Curve::generator());
    }

    #[test]
    fn msm_23_is_6_over_field_elements() {
        let c: [FrElement; 1] = [FrElement::from(3)];
        let hiding = [G1::new(2)];
        assert_eq!(msm(&c, &hiding), G1::new(6));
    }

    #[test]
    fn msm_23_is_6_over_elliptic_curves() {
        let c: [FrElement; 1] = [FrElement::from(3)];
        let g = BLS12381Curve::generator();
        let hiding = [g.operate_with_self(2_u16)];
        assert_eq!(msm(&c, &hiding), g.operate_with_self(6_u16));
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18_over_field_elements() {
        let c: [FrElement; 2] = [FrElement::from(2), FrElement::from(3)];
        let hiding = [G1::new(3), G1::new(4)];
        assert_eq!(msm(&c, &hiding), G1::new(18));
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18_over_elliptic_curves() {
        let c: [FrElement; 2] = [FrElement::from(2), FrElement::from(3)];
        let g = BLS12381Curve::generator();
        let hiding = [g.operate_with_self(3_u16), g.operate_with_self(4_u16)];
        assert_eq!(msm(&c, &hiding), g.operate_with_self(18_u16));
    }

    #[test]
    fn msm_with_empty_input_over_field_elements() {
        let c: [FrElement; 0] = [];
        let hiding: [G1; 0] = [];
        assert_eq!(msm(&c, &hiding), G1::new(0));
    }

    #[test]
    fn msm_with_empty_c_is_none_over_elliptic_curves() {
        let c: [FrElement; 0] = [];
        let hiding: [G1; 0] = [];
        assert_eq!(
            msm(&c, &hiding),
            ShortWeierstrassProjectivePoint::neutral_element()
        );
    }
}
