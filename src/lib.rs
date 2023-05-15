#![allow(unused_variables)]
#![allow(clippy::not_unsafe_ptr_arg_deref)]

pub mod commitments;
pub mod compress;
pub mod math;
pub mod sqrt;
pub mod srs;
pub mod utils;

use crate::compress::{compress_g1_point, decompress_g1_point};
use crate::math::cyclic_group::IsGroup;
pub use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::default_types::{
    FrConfig, FrElement, FrField, MODULUS,
};
use crate::math::elliptic_curve::short_weierstrass::curves::bls12_381::field_extension::LevelOneResidue;
use crate::math::elliptic_curve::traits::FromAffine;
use crate::math::elliptic_curve::traits::IsEllipticCurve;
use crate::math::field::extensions::quadratic::QuadraticExtensionField;
use commitments::{
    kzg::KateZaveruchaGoldberg, kzg::StructuredReferenceString, traits::IsCommitmentScheme,
};
use core::ptr::null_mut;
use math::polynomial::Polynomial;
use math::{
    elliptic_curve::short_weierstrass::{
        curves::bls12_381::{
            curve::BLS12381Curve, field_extension::BLS12381PrimeField, pairing::BLS12381AtePairing,
            twist::BLS12381TwistCurve,
        },
        point::ShortWeierstrassProjectivePoint,
    },
    field::element::FieldElement,
    msm::g1_lincomb,
    traits::ByteConversion,
};
use std::marker;

pub type G1 = ShortWeierstrassProjectivePoint<BLS12381Curve>;
pub type G1Point = ShortWeierstrassProjectivePoint<BLS12381Curve>;
pub type G2Point = <BLS12381TwistCurve as IsEllipticCurve>::PointRepresentation;
pub type KZG = KateZaveruchaGoldberg<FrField, BLS12381AtePairing>;
pub type BLS12381FieldElement = FieldElement<BLS12381PrimeField>;
pub type BLS12381TwistCurveFieldElement = FieldElement<QuadraticExtensionField<LevelOneResidue>>;
#[allow(clippy::upper_case_acronyms)]
type FE = FrElement;

#[derive(Debug, PartialEq)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub enum C_KZG_RET {
    /// Success
    C_KZG_OK = 0,
    /// The supplied data is invalid in some way
    C_KZG_BADARGS,
    /// Internal error - this should never occur
    C_KZG_ERROR,
    /// Could not allocate memory
    C_KZG_MALLOC,
}

/** The domain separator for the Fiat-Shamir protocol. */
pub const FIAT_SHAMIR_PROTOCOL_DOMAIN: [u8; 16] = *b"FSBLOBVERIFY_V1_";

pub const RANDOM_CHALLENGE_KZG_BATCH_DOMAIN: [u8; 16] = *b"RCKZGBATCH___V1_";

/** Length of the domain strings above. */
pub const DOMAIN_STR_LENGTH: usize = 16;

/** The number of bytes in a KZG commitment. */
pub const BYTES_PER_COMMITMENT: usize = 48;

/** The number of bytes in a KZG proof. */
pub const BYTES_PER_PROOF: usize = 48;

/** The number of bytes in a BLS scalar field element. */
pub const BYTES_PER_FIELD_ELEMENT: usize = 32;

pub const FIELD_ELEMENTS_PER_BLOB: usize = 4096;

pub const TRUSTED_SETUP_NUM_G1_POINTS: usize = FIELD_ELEMENTS_PER_BLOB;

/** The number of bytes in a blob. */
pub const BYTES_PER_BLOB: usize = FIELD_ELEMENTS_PER_BLOB * BYTES_PER_FIELD_ELEMENT;

pub const CHALLENGE_INPUT_SIZE: usize =
    DOMAIN_STR_LENGTH + 16 + BYTES_PER_BLOB + BYTES_PER_COMMITMENT;

pub const BYTES_PER_G1_POINT: usize = 48;
pub const BYTES_PER_G2_POINT: usize = 96;

pub const NUM_G1_POINTS: usize = FIELD_ELEMENTS_PER_BLOB;
/// Number of G2 points required for the kzg trusted setup.
/// 65 is fixed and is used for providing multiproofs up to 64 field elements.
pub const NUM_G2_POINTS: usize = 65;

pub type Bytes32 = [u8; 32];
pub type Bytes48 = [u8; 48];
pub type KZGCommitment = [u8; 48];
pub type KZGProof = [u8; 48];
pub type Blob = [u8; BYTES_PER_BLOB];

#[allow(non_camel_case_types)]
pub type limb_t = u64;

#[derive(Default)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_fr {
    l: [limb_t; 256 / 8 / core::mem::size_of::<limb_t>()],
}

#[derive(Default, Copy, Clone)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_fp {
    l: [limb_t; 384 / 8 / core::mem::size_of::<limb_t>()],
}

#[derive(Default, Copy, Clone)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p1 {
    pub x: blst_fp,
    pub y: blst_fp,
    pub z: blst_fp,
}

#[derive(Default)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p1_affine {
    pub x: blst_fp,
    pub y: blst_fp,
}

/* 0 is "real" part, 1 is "imaginary" */
#[derive(Default, Copy, Clone)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_fp2 {
    fp: [blst_fp; 2],
}

#[allow(non_camel_case_types)]
pub type g1_t = blst_p1;
/**< Internal G1 group element type. */
#[allow(non_camel_case_types)]
pub type g2_t = blst_p2;
/**< Internal G2 group element type. */
#[allow(non_camel_case_types)]
pub type fr_t = blst_fr;
/**< Internal Fr field element type. */

#[derive(Default, Clone, Copy)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p2 {
    pub x: blst_fp2,
    pub y: blst_fp2,
    pub z: blst_fp2,
}

#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p2_affine {
    pub x: blst_fp2,
    pub y: blst_fp2,
}

//typedef struct { limb_t l[256/8/sizeof(limb_t)]; } blst_fr;
//typedef struct { limb_t l[384/8/sizeof(limb_t)]; } blst_fp;
//typedef struct { blst_fp x, y, z; } blst_p1;
//typedef struct { blst_fp x, y; } blst_p1_affine;

#[allow(non_camel_case_types)]
#[repr(C)]
/// Stores the setup and parameters needed for performing FFTs.
pub struct FFTSettings<'a> {
    /** The maximum size of FFT these settings support, a power of 2. */
    pub max_width: u64,
    /**
     * We fix a given primitive roots of unity w of order `max_width` via
     * `SCALE2_ROOT_OF_UNITY`. Then `expanded_roots_of_unity[i]` == w^i and
     * `reverse_roots_of_unity[i]` == w^{-i}. Unusually, both
     * `expanded_roots_of_unity` and `reverse_roots_of_unity` have length
     * `max_width + 1`. By the above, `expanded_roots_of_unity[max_width] ==
     * expanded_roots_of_unity[0] == 1` and similarly for
     * `reverse_roots_of_unity`. The redundant element is just there to simplify
     * index calculations in some formulas.
     */

    /** Ascending powers of the root of unity, length `max_width + 1`. */
    pub expanded_roots_of_unity: *mut fr_t,
    /** Descending powers of the root of unity, length `max_width + 1`. */
    pub reverse_roots_of_unity: *mut fr_t,
    /** Powers of the root of unity in bit-reversal permutation order, length
     * `max_width`. */
    pub roots_of_unity: *mut fr_t,

    _marker: marker::PhantomData<&'a *mut fr_t>,
}

impl<'a> Default for FFTSettings<'a> {
    fn default() -> Self {
        Self {
            max_width: 0,
            expanded_roots_of_unity: null_mut(),
            reverse_roots_of_unity: null_mut(),
            roots_of_unity: null_mut(),
            _marker: marker::PhantomData,
        }
    }
}

#[derive(Clone)]
#[allow(non_camel_case_types)]
#[repr(C)]
pub struct KZGSettings<'a> {
    /** The corresponding settings for performing FFTs. */
    pub fs: *mut FFTSettings<'a>,
    /** G1 group elements from the trusted setup,
     * in Lagrange form bit-reversal permutation. */
    pub g1_values: *mut g1_t,
    /** G2 group elements from the trusted setup;
     * both arrays have `FIELD_ELEMENTS_PER_BLOB` elements. */
    pub g2_values: *mut g2_t,

    _marker: marker::PhantomData<&'a *mut FFTSettings<'a>>,
    _marker2: marker::PhantomData<&'a *mut g1_t>,
    _marker3: marker::PhantomData<&'a *mut g2_t>,
}

impl<'a> Default for KZGSettings<'a> {
    fn default() -> Self {
        Self {
            fs: null_mut(),
            g1_values: null_mut(),
            g2_values: null_mut(),
            _marker: marker::PhantomData,
            _marker2: marker::PhantomData,
            _marker3: marker::PhantomData,
        }
    }
}

#[no_mangle]
/// Convert a blob to a KZG commitment.
///
/// # Params
///
/// * `out` -  The resulting commitment
/// * `blob` - The blob representing the polynomial to be committed to
/// * `s`    - The trusted setup
pub extern "C" fn blob_to_kzg_commitment(
    out: *mut KZGCommitment,
    blob: *const Blob,
    s: *const KZGSettings,
) -> C_KZG_RET {
    let s_struct = unsafe { (*s).clone() };
    let input_blob: [u8; BYTES_PER_BLOB] =
        unsafe { std::slice::from_raw_parts(blob, BYTES_PER_BLOB)[0] };

    let Ok(polynomial) = utils::blob_to_polynomial(&input_blob) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let srs = kzgsettings_to_structured_reference_string(&s_struct);
    let kzg = KZG::new(srs);
    let commitment: ShortWeierstrassProjectivePoint<BLS12381Curve> = kzg.commit(&polynomial);
    let Ok(commitment_bytes) = compress_g1_point(&commitment) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    unsafe {
        std::ptr::copy(
            commitment_bytes.as_ptr(),
            out as *mut u8,
            BYTES_PER_COMMITMENT,
        );
    }
    C_KZG_RET::C_KZG_OK
}

///
/// Compute KZG proof for polynomial in Lagrange form at position z.
///
/// # Params
///
/// `proof_out` - The combined proof as a single G1 element
/// `y_out` - The evaluation of the polynomial at the evaluation point z
/// `blob` - The blob (polynomial) to generate a proof for
/// `z` - The generator z-value for the evaluation points
/// `s` - The trusted setup
///
/// # Return
///
/// Value of type `C_KZG_RET` indicating error status.
#[no_mangle]
pub extern "C" fn compute_kzg_proof(
    proof_out: *mut KZGProof,
    y_out: *mut Bytes32,
    blob: *const Blob,
    z_bytes: *const Bytes32,
    s: *const KZGSettings,
) -> C_KZG_RET {
    let z_slice = unsafe { *z_bytes };
    let s_struct = unsafe { (*s).clone() };
    let input_blob: [u8; BYTES_PER_BLOB] =
        unsafe { std::slice::from_raw_parts(blob, BYTES_PER_BLOB)[0] };

    let Ok(polynomial) = utils::blob_to_polynomial(&input_blob) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(fr_z) = FrElement::from_bytes_be(&z_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let fr_y: FE = polynomial.evaluate(&fr_z);
    let y_out_slice: [u8; 32] = fr_y.to_bytes_be().try_into().unwrap();

    let srs = kzgsettings_to_structured_reference_string(&s_struct);
    let kzg = KZG::new(srs);
    let proof = kzg.open(&fr_z, &fr_y, &polynomial);
    let Ok(compressed_proof) = compress_g1_point(&proof) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    unsafe {
        std::ptr::copy(
            compressed_proof.as_ptr(),
            proof_out as *mut u8,
            BYTES_PER_PROOF,
        );
        std::ptr::copy(y_out_slice.as_ptr(), y_out as *mut u8, 32);
    }

    C_KZG_RET::C_KZG_OK
}

/// Given a blob and a commitment, return the KZG proof that is used to verify
/// it against the commitment. This function does not verify that the commitment
/// is correct with respect to the blob.
///
/// # Params
///
/// `out` - The resulting proof
/// `blob` - A blob representing a polynomial
/// `commitment_bytes` - Commitment to verify
/// `s` - The trusted setup
///
/// # Return
///
/// Value of type `C_KZG_RET` indicating error status.
#[no_mangle]
pub extern "C" fn compute_blob_kzg_proof(
    out: *mut KZGProof,
    blob: *const Blob,
    commitment_bytes: *const Bytes48,
    s: *const KZGSettings,
) -> C_KZG_RET {
    let mut commitment_slice = unsafe { *commitment_bytes };
    let input_blob: [u8; BYTES_PER_BLOB] =
        unsafe { std::slice::from_raw_parts(blob, BYTES_PER_BLOB)[0] };
    let s_struct = unsafe { (*s).clone() };

    // Do conversions first to fail fast, compute_challenge is expensive
    let Ok(commitment_g1) = decompress_g1_point(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };
    let Ok(polynomial) = utils::blob_to_polynomial(&input_blob) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    // get the point fr_z (i.e. "x") from fr where evaluate the polynomial
    // Compute the challenge for the given blob/commitment
    let Ok(fr_z) = utils::compute_challenge(
        &input_blob,
        &commitment_g1,
    ) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let fr_y: FE = polynomial.evaluate(&fr_z);
    let srs = kzgsettings_to_structured_reference_string(&s_struct);
    let kzg = KZG::new(srs);
    let proof = kzg.open(&fr_z, &fr_y, &polynomial);
    let Ok(compressed_proof) = compress_g1_point(&proof) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    unsafe {
        std::ptr::copy(compressed_proof.as_ptr(), out as *mut u8, BYTES_PER_PROOF);
    }

    C_KZG_RET::C_KZG_OK
}

#[no_mangle]
pub extern "C" fn verify_kzg_proof(
    ok: *mut bool,
    commitment_bytes: *const Bytes48,
    z_bytes: *const Bytes32,
    y_bytes: *const Bytes32,
    proof_bytes: *const Bytes48,
    s: *const KZGSettings,
) -> C_KZG_RET {
    unsafe {
        *ok = false;
    }
    let mut commitment_slice = unsafe { *commitment_bytes };
    let z_slice = unsafe { *z_bytes };
    let y_slice = unsafe { *y_bytes };
    let mut proof_slice = unsafe { *proof_bytes };
    let s_struct = unsafe { (*s).clone() };

    let Ok(commitment_g1) = decompress_g1_point(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(z_fr) = FrElement::from_bytes_be(&z_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(y_fr) = FrElement::from_bytes_be(&y_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(proof_g1) = decompress_g1_point(&mut proof_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    // FIXME: We should not use create_src() for this instantiation.
    let kzg = KZG::new(utils::create_srs());
    let ret = kzg.verify(&z_fr, &y_fr, &commitment_g1, &proof_g1);

    if ret {
        unsafe {
            *ok = true;
        }
    }

    C_KZG_RET::C_KZG_OK
}

#[no_mangle]
pub extern "C" fn verify_blob_kzg_proof(
    ok: *mut bool,
    blob: *const Blob,
    commitment_bytes: *const Bytes48,
    proof_bytes: *const Bytes48,
    s: *const KZGSettings,
) -> C_KZG_RET {
    unsafe {
        *ok = false;
    }
    let input_blob: [u8; BYTES_PER_BLOB] =
        unsafe { std::slice::from_raw_parts(blob, BYTES_PER_BLOB)[0] };

    let mut commitment_slice = unsafe { *commitment_bytes };
    let mut proof_slice = unsafe { *proof_bytes };
    let s_struct = unsafe { (*s).clone() };

    let Ok(commitment_g1) = decompress_g1_point(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };
    let Ok(proof_g1) = decompress_g1_point(&mut proof_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };
    let Ok(polynomial) = utils::blob_to_polynomial(&input_blob) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(evaluation_challenge_fr) = utils::compute_challenge(
        &input_blob,
        &commitment_g1,
    ) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let y_fr: FE = polynomial.evaluate(&evaluation_challenge_fr);

    let srs = kzgsettings_to_structured_reference_string(&s_struct);
    let kzg = KZG::new(srs);
    let ret = kzg.verify(&evaluation_challenge_fr, &y_fr, &commitment_g1, &proof_g1);

    if ret {
        unsafe {
            *ok = true;
        }
    }

    C_KZG_RET::C_KZG_OK
}

/// Given a list of blobs and blob KZG proofs, verify that they correspond to the
/// provided commitments.
///
/// This function assumes that `n` is trusted and that all input arrays
/// contain `n` elements. `n` should be the actual size of the arrays and not
/// read off a length field in the protocol.
///
/// This function accepts if called with `n==0`.
///
/// # Params
///
/// * `ok` - True if the proofs are valid, otherwise false
/// * `blobs` - Array of blobs to verify
/// * `commitments_bytes` - Array of commitments to verify
/// * `proofs_bytes` - Array of proofs used for verification
/// * `n` - The number of blobs/commitments/proofs
/// * `s` - The trusted setup
#[no_mangle]
pub extern "C" fn verify_blob_kzg_proof_batch(
    ok: *mut bool,
    blobs: *const Blob,
    commitments_bytes: *const Bytes48,
    proofs_bytes: *const Bytes48,
    n: usize,
    s: *const KZGSettings,
) -> C_KZG_RET {
    unsafe {
        *ok = false;
    }

    match n {
        0 => {
            unsafe {
                *ok = false;
            }
            C_KZG_RET::C_KZG_OK
        }
        1 => verify_blob_kzg_proof(ok, blobs, commitments_bytes, proofs_bytes, s),
        n => {
            let s_struct = unsafe { (*s).clone() };
            let commitments_slice_array =
                unsafe { std::slice::from_raw_parts(commitments_bytes, n) };
            let mut commitments_g1_vec = Vec::<G1Point>::new();

            // arrays of blobs
            let blobs_slice_array = unsafe { std::slice::from_raw_parts(blobs, n) };
            let mut polynomial_vec = Vec::<Polynomial<FrElement>>::new();

            // proofs
            let proofs_slice_array = unsafe { std::slice::from_raw_parts(proofs_bytes, n) };

            let mut challenges_vec = Vec::<FE>::new();
            let mut evaluations_vec = Vec::<FE>::new();
            let mut proofs_vec = Vec::<G1Point>::new();

            for i in 0..n {
                // fill the vector of commitments
                let mut g1_point_compressed = commitments_slice_array[i];
                let Ok(commitment_g1) = decompress_g1_point(&mut g1_point_compressed) else {
                    return C_KZG_RET::C_KZG_ERROR;
                };

                // convert blobs_slice_array to polynomials
                let Ok(polynomial) = utils::blob_to_polynomial(&blobs_slice_array[i]) else {
                    return C_KZG_RET::C_KZG_ERROR;
                };

                // get list of challenges
                // Compute the challenge for the given blob/commitment
                let Ok(fr_z) = utils::compute_challenge(
                    &blobs_slice_array[i],
                    &commitment_g1,
                ) else {
                    return C_KZG_RET::C_KZG_ERROR;
                };

                // get list of evaluations
                let fr_y: FE = polynomial.evaluate(&fr_z);
                evaluations_vec.push(fr_y);

                challenges_vec.push(fr_z);
                commitments_g1_vec.push(commitment_g1);
                polynomial_vec.push(polynomial);

                let mut proof_slice = proofs_slice_array[i];
                let Ok(proof_g1) = decompress_g1_point(&mut proof_slice) else {
                    return C_KZG_RET::C_KZG_ERROR;
                };
                proofs_vec.push(proof_g1);
            }
            let Ok(ret_verify) = verify_kzg_proof_batch(
                &commitments_g1_vec,
                &challenges_vec,
                &evaluations_vec,
                &proofs_vec,
                n,
                &s_struct,
            ) else {
                return C_KZG_RET::C_KZG_ERROR;
            };
            unsafe {
                *ok = ret_verify;
            }

            C_KZG_RET::C_KZG_OK
        }
    }
}

///
/// Helper function for verify_blob_kzg_proof_batch(): actually perform the
/// verification.
///
/// Remark: This function assumes that `n` is trusted and that all input arrays
///     contain `n` elements. `n` should be the actual size of the arrays and not
///     read off a length field in the protocol.
///
/// Remark: This function only works for `n > 0`.
///
/// # Params
///
/// * `ok` - True if the proofs are valid, otherwise false
/// * `commitments_g1` - Array of commitments to verify
/// * `zs_fr` - Array of evaluation points for the KZG proofs
/// * `ys_fr` - Array of evaluation results for the KZG proofs
/// * `proofs_g1` - Array of proofs used for verification
/// * `n` - The number of blobs/commitments/proofs
/// * `s` - The trusted setup
///
/// # Returns
///
/// * `Ok(true)` if the proofs are valid, otherwise `Ok(false)`
fn verify_kzg_proof_batch(
    commitments_g1: &[G1Point],
    zs_fr: &[FE],
    ys_fr: &[FE],
    proofs_g1: &[G1Point],
    n: usize,
    s: &KZGSettings,
) -> Result<bool, Vec<u8>> {
    let r_powers: Vec<_> = utils::compute_r_powers(commitments_g1, zs_fr, ys_fr, proofs_g1)?
        .iter()
        .map(|x| x.representative())
        .collect();

    let mut c_minus_y = Vec::new();
    let mut r_times_z = Vec::new();
    // Compute \sum r^i * Proof_i
    let proof_lincomb = g1_lincomb(proofs_g1, &r_powers);
    let g = BLS12381Curve::generator();

    for i in 0..n {
        let ys_encrypted = g.operate_with_self(ys_fr[i].representative());
        // Get C_i - [y_i]
        let c_i = commitments_g1[i].clone();
        let c_minus_y_elem = c_i.operate_with(&ys_encrypted.neg());
        c_minus_y.push(c_minus_y_elem);

        // Get r^i * z_i
        let r_i = FE::from(&r_powers[i]);
        let zs_fr_i = zs_fr[i].clone();
        let r_times_z_elem = r_i * zs_fr_i;
        r_times_z.push(r_times_z_elem.representative());
    }

    //  Get \sum r^i z_i Proof_i
    let proof_z_lincomb = g1_lincomb(proofs_g1, &r_times_z);
    //
    // Get \sum r^i (C_i - [y_i])
    let c_minus_y_lincomb = g1_lincomb(&c_minus_y, &r_powers);

    // Get c_minus_y_lincomb + proof_z_lincomb
    let rhs_g1 = c_minus_y_lincomb.operate_with(&proof_z_lincomb);

    let srs = kzgsettings_to_structured_reference_string(s);
    let kzg = KZG::new(srs);
    Ok(kzg.verify(&FE::zero(), &FE::zero(), &rhs_g1, &proof_z_lincomb))
}

/* TODO: implement
C_KZG_RET load_trusted_setup_file(KZGSettings *out, FILE *in);
C_KZG_RET load_trusted_setup(
    KZGSettings *out,
    const uint8_t *g1_bytes, /* n1 * 48 bytes */
    size_t n1,
    const uint8_t *g2_bytes, /* n2 * 96 bytes */
    size_t n2
);
*/

pub fn vecs_to_structured_reference_string(
    g1_points: &[G1],
    g2_points: &[G2Point],
) -> StructuredReferenceString<G1, G2Point> {
    StructuredReferenceString::<G1, G2Point>::new(g1_points, g2_points)
}

pub fn kzgsettings_to_structured_reference_string(
    s: &KZGSettings,
) -> StructuredReferenceString<G1, G2Point> {
    let g1_values = s.g1_values;
    let g2_values = s.g2_values;

    let g1_points_slice: [blst_p1; NUM_G1_POINTS] = unsafe { *g1_values.cast() };
    let g2_points_slice: [blst_p2; NUM_G2_POINTS] = unsafe { *g2_values.cast() };

    let g1_points: Vec<G1> = g1_points_slice
        .iter()
        .map(|point| {
            let x_field = point
                .x
                .l
                .iter()
                .flat_map(|e| e.to_be_bytes())
                .collect::<Vec<u8>>();
            let x = BLS12381FieldElement::from_bytes_be(&x_field).unwrap();

            let y_field = point
                .y
                .l
                .iter()
                .flat_map(|e| e.to_be_bytes())
                .collect::<Vec<u8>>();
            let y = BLS12381FieldElement::from_bytes_be(&y_field).unwrap();
            G1::from_affine(x, y).unwrap()
        })
        .collect();

    let g2_points: Vec<G2Point> = g2_points_slice
        .iter()
        .map(|point| {
            let [x0, x1] = point.x.fp;
            let [y0, y1] = point.y.fp;
            //let z = point.z;

            let x0_field = BLS12381FieldElement::from_bytes_be(
                &x0.l
                    .iter()
                    .flat_map(|e| e.to_be_bytes())
                    .collect::<Vec<u8>>(),
            )
            .unwrap();
            let x1_field = BLS12381FieldElement::from_bytes_be(
                &x1.l
                    .iter()
                    .flat_map(|e| e.to_be_bytes())
                    .collect::<Vec<u8>>(),
            )
            .unwrap();
            let y0_field = BLS12381FieldElement::from_bytes_be(
                &y0.l
                    .iter()
                    .flat_map(|e| e.to_be_bytes())
                    .collect::<Vec<u8>>(),
            )
            .unwrap();
            let y1_field = BLS12381FieldElement::from_bytes_be(
                &y1.l
                    .iter()
                    .flat_map(|e| e.to_be_bytes())
                    .collect::<Vec<u8>>(),
            )
            .unwrap();

            let x = BLS12381TwistCurveFieldElement::new([x0_field, x1_field]);
            let y = BLS12381TwistCurveFieldElement::new([y0_field, y1_field]);
            G2Point::from_affine(x, y).unwrap()
        })
        .collect();

    StructuredReferenceString::<G1, G2Point>::new(&g1_points, &g2_points)
}

#[cfg(test)]
mod tests {
    use crate::commitments::traits::IsCommitmentScheme;
    use crate::compress::{compress_g1_point, decompress_g1_point};
    use crate::math::{
        cyclic_group::IsGroup, field::element::FieldElement, polynomial::Polynomial,
        traits::ByteConversion,
    };
    use crate::srs::load_trusted_setup_file_to_g1_points_and_g2_points;
    use crate::utils::polynomial_to_blob_with_size;
    use crate::{
        blst_fr, blst_p1, blst_p2, compute_kzg_proof, fr_t, verify_blob_kzg_proof_batch,
        verify_kzg_proof, Blob, Bytes32, Bytes48, FFTSettings, FrElement, G1Point, KZGProof,
        KZGSettings, BYTES_PER_BLOB, C_KZG_RET, FE,
    };

    #[test]
    fn test_compute_kzg_proof() {
        // Test this case:
        // polinomial: 1 cte
        // evaluation fr: 1
        // Expected result:
        // proof: inf
        // y_out: 1

        // output buffers
        let mut proof_out: KZGProof = [0; 48];
        let mut y_out: Bytes32 = [0; 32];

        // assign poly 1
        let polynomial = Polynomial::<FrElement>::new(&[FieldElement::one()]);
        let blob = polynomial_to_blob_with_size(&polynomial).unwrap();
        // z = 1
        let z_bytes: Bytes32 = FE::from(1).to_bytes_be().try_into().unwrap();

        let mut zero_fr_expanded_roots_of_unity = blst_fr::default();
        let mut zero_fr_reverse_roots_of_unity = blst_fr::default();
        let mut zero_fr_roots_of_unity = blst_fr::default();

        let mut fft_settings = FFTSettings {
            max_width: 8,
            expanded_roots_of_unity: &mut zero_fr_expanded_roots_of_unity as *mut fr_t,
            reverse_roots_of_unity: &mut zero_fr_reverse_roots_of_unity as *mut fr_t,
            roots_of_unity: &mut zero_fr_roots_of_unity as *mut fr_t,
            _marker: std::marker::PhantomData,
        };

        let mut zer_fr_g1_values = blst_p1::default();
        let mut zer_fr_g2_values = blst_p2::default();

        let s = KZGSettings {
            fs: &mut fft_settings as *mut FFTSettings,
            g1_values: &mut zer_fr_g1_values as *mut blst_p1,
            g2_values: &mut zer_fr_g2_values as *mut blst_p2,
            _marker: std::marker::PhantomData,
            _marker2: std::marker::PhantomData,
            _marker3: std::marker::PhantomData,
        };

        let ret = compute_kzg_proof(
            &mut proof_out as *mut KZGProof,
            &mut y_out as *mut Bytes32,
            &blob as *const Blob,
            &z_bytes as *const Bytes32,
            &s as *const KZGSettings,
        );

        // assert ret function
        let ok_enum_kzg = C_KZG_RET::C_KZG_OK;
        assert_eq!(ret, ok_enum_kzg);

        // y evaluation is one
        let one_fr = FE::one();
        let y_out_fr = FE::from_bytes_be(&y_out).unwrap();
        assert_eq!(one_fr, y_out_fr);

        // proof is inf
        let p = decompress_g1_point(&mut proof_out).unwrap();
        let inf = G1Point::neutral_element();
        assert_eq!(p, inf);

        let kzg = crate::KZG::new(crate::utils::create_srs());
        let commitment = kzg.commit(&polynomial);
        let commitment_bytes = compress_g1_point(&commitment).unwrap();

        let mut ok = false;
        let ret_verify = verify_kzg_proof(
            &mut ok as *mut bool,
            &commitment_bytes as *const Bytes48,
            &z_bytes as *const Bytes32,
            &y_out as *const Bytes32,
            &proof_out as *const Bytes48,
            &s as *const KZGSettings,
        );

        assert_eq!(ret_verify, ok_enum_kzg);

        // FIXME: make blobs useful
        let blobs: Blob = [0; BYTES_PER_BLOB];

        // verify blob as a batch
        ok = false;
        verify_blob_kzg_proof_batch(
            &mut ok as *mut bool,
            &blobs as *const Blob,
            &commitment_bytes as *const Bytes48,
            &proof_out as *const Bytes48,
            1,
            &s as *const KZGSettings,
        );
    }

    #[test]
    fn test_read_srs() {
        let (g1_points, g2_points) =
            load_trusted_setup_file_to_g1_points_and_g2_points("test/trusted_setup_4.txt").unwrap();

        let settings = super::vecs_to_structured_reference_string(&g1_points, &g2_points);
    }
}
