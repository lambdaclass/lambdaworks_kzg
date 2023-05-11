use crate::math::cyclic_group::IsGroup;
use crate::math::unsigned_integer::element::UnsignedInteger;
use crate::math::unsigned_integer::traits::IsUnsignedInteger;
use crate::G1Point;

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

/// Calculate a linear combination of G1 group elements.
///
/// Calculates `[coeffs_0]p_0 + [coeffs_1]p_1 + ... + [coeffs_n]p_n`
/// where `n` is `len - 1`.
///
/// This function computes the result using Pippenger's algorithm.
pub fn g1_lincomb(p: &[G1Point], coeff: &[UnsignedInteger<4>]) -> G1Point {
    msm_pip(coeff, p)
}

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
pub fn msm_pip<const NUM_LIMBS: usize, G>(cs: &[UnsignedInteger<NUM_LIMBS>], hidings: &[G]) -> G
where
    G: IsGroup,
{
    debug_assert_eq!(
        cs.len(),
        hidings.len(),
        "Slices `cs` and `hidings` must be of the same length to compute `msm`."
    );
    // When input is small enough, windows of length 2 seem faster than 1.
    const MIN_WINDOWS: usize = 2;
    const SCALE_FACTORS: (usize, usize) = (4, 5);

    // We approximate the optimum window size with: f(n) = k * log2(n), where k is a scaling factor
    let window_size =
        ((usize::BITS - cs.len().leading_zeros() - 1) as usize * SCALE_FACTORS.0) / SCALE_FACTORS.1;
    msm_with(cs, hidings, MIN_WINDOWS.min(window_size))
}

pub fn msm_with<const NUM_LIMBS: usize, G>(
    cs: &[UnsignedInteger<NUM_LIMBS>],
    hidings: &[G],
    window_size: usize,
) -> G
where
    G: IsGroup,
{
    debug_assert!(window_size < usize::BITS as usize);
    // The number of windows of size `s` is ceil(lambda/s).
    let num_windows = (64 * NUM_LIMBS - 1) / window_size + 1;

    // We define `buckets` outside of the loop so we only have to allocate once, and reuse it.
    //
    // This line forces a heap allocation which might be undesired. We can define this buckets
    // variable in the Pippenger struct to only allocate once, but use a bit of extra memory.
    // If we accept a const window_size, we could make it an array instaed of a vector
    // avoiding the heap allocation. We should be aware if that might be too agressive for
    // the stack and cause a potential stack overflow.
    let mut buckets = vec![G::neutral_element(); (1 << window_size) - 1];

    (0..num_windows)
        .rev()
        .map(|window_idx| {
            // Put in the right bucket the corresponding ps[i] for the current window.
            cs.iter().zip(hidings).for_each(|(k, p)| {
                // We truncate the number to the least significative limb.
                // This is ok because window_size < usize::BITS.
                let window_unmasked = (k >> (window_idx * window_size)).limbs[NUM_LIMBS - 1];
                let m_ij = window_unmasked & ((1 << window_size) - 1);
                if m_ij != 0 {
                    let idx = (m_ij - 1) as usize;
                    buckets[idx] = buckets[idx].operate_with(p);
                }
            });

            // Do the reduction step for the buckets.
            buckets
                .iter_mut()
                .rev()
                .scan(G::neutral_element(), |m, b| {
                    *m = m.operate_with(b); // Reduction step.
                    *b = G::neutral_element(); // Cleanup bucket slot to reuse in the next window.
                    Some(m.clone())
                })
                .reduce(|g, m| g.operate_with(&m))
                .unwrap_or_else(G::neutral_element)
        })
        .reduce(|t, g| t.operate_with_self(1_u64 << window_size).operate_with(&g))
        .unwrap_or_else(G::neutral_element)
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
