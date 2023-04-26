use crate::math::cyclic_group::IsGroup;
use crate::math::unsigned_integer::traits::IsUnsignedInteger;

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
pub fn msm<C: IsUnsignedInteger, T>(cs: &[C], hidings: &[T]) -> T
where
    T: IsGroup,
{
    debug_assert_eq!(
        cs.len(),
        hidings.len(),
        "Slices `cs` and `hidings` must be of the same length to compute `msm`."
    );
    cs.iter()
        .zip(hidings.iter())
        .map(|(&c, h)| h.operate_with_self(c))
        .reduce(|acc, x| acc.operate_with(&x))
        .unwrap_or_else(T::neutral_element)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::elliptic_curve::short_weierstrass::curves::test_curve_1::TestCurve1;
    use crate::math::elliptic_curve::short_weierstrass::point::ShortWeierstrassProjectivePoint;
    use crate::math::elliptic_curve::traits::IsEllipticCurve;
    use crate::math::field::fields::u64_prime_field::U64FieldElement;

    const ORDER_R: u64 = 5;
    type FE = U64FieldElement<ORDER_R>;

    #[test]
    fn msm_11_is_1_over_elliptic_curves() {
        let c: [u64; 1] = [1];
        let hiding = [TestCurve1::generator()];
        assert_eq!(msm(&c, &hiding), TestCurve1::generator());
    }

    #[test]
    fn msm_23_is_6_over_field_elements() {
        let c: [u64; 1] = [3];
        let hiding = [FE::new(2)];
        assert_eq!(msm(&c, &hiding), FE::new(6));
    }

    #[test]
    fn msm_23_is_6_over_elliptic_curves() {
        let c: [u64; 1] = [3];
        let g = TestCurve1::generator();
        let hiding = [g.operate_with_self(2_u16)];
        assert_eq!(msm(&c, &hiding), g.operate_with_self(6_u16));
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18_over_field_elements() {
        let c: [u64; 2] = [2, 3];
        let hiding = [FE::new(3), FE::new(4)];
        assert_eq!(msm(&c, &hiding), FE::new(18));
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18_over_elliptic_curves() {
        let c: [u64; 2] = [2, 3];
        let g = TestCurve1::generator();
        let hiding = [g.operate_with_self(3_u16), g.operate_with_self(4_u16)];
        assert_eq!(msm(&c, &hiding), g.operate_with_self(18_u16));
    }

    #[test]
    fn msm_with_empty_input_over_field_elements() {
        let c: [u64; 0] = [];
        let hiding: [FE; 0] = [];
        assert_eq!(msm(&c, &hiding), FE::new(0));
    }

    #[test]
    fn msm_with_empty_c_is_none_over_elliptic_curves() {
        let c: [u64; 0] = [];
        let hiding: [ShortWeierstrassProjectivePoint<TestCurve1>; 0] = [];
        assert_eq!(
            msm(&c, &hiding),
            ShortWeierstrassProjectivePoint::neutral_element()
        );
    }
}
