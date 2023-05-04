#![allow(unused_variables)]
#![allow(clippy::not_unsafe_ptr_arg_deref)]

pub mod commitments;
pub mod math;
pub mod utils;

use commitments::{
    kzg::{FrElement, FrField, KateZaveruchaGoldberg},
    traits::IsCommitmentScheme,
};
use math::{
    elliptic_curve::short_weierstrass::{
        curves::bls12_381::{
            curve::BLS12381Curve, field_extension::BLS12381PrimeField, pairing::BLS12381AtePairing,
            twist::BLS12381TwistCurve,
        },
        point::ShortWeierstrassProjectivePoint,
    },
    field::element::FieldElement,
    traits::ByteConversion,
};
use std::marker;

use crate::utils::compress_g1_point;

pub type G1Point = ShortWeierstrassProjectivePoint<BLS12381Curve>;
pub type G2Point = ShortWeierstrassProjectivePoint<BLS12381TwistCurve>;
pub type KZG = KateZaveruchaGoldberg<FrField, BLS12381AtePairing>;

pub type BLS12381FieldElement = FieldElement<BLS12381PrimeField>;

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

/** The number of bytes in a KZG commitment. */
pub const BYTES_PER_COMMITMENT: usize = 48;

/** The number of bytes in a KZG proof. */
pub const BYTES_PER_PROOF: usize = 48;

/** The number of bytes in a BLS scalar field element. */
pub const BYTES_PER_FIELD_ELEMENT: usize = 32;

pub const FIELD_ELEMENTS_PER_BLOB: usize = 4096;

/** The number of bytes in a blob. */
pub const BYTES_PER_BLOB: usize = FIELD_ELEMENTS_PER_BLOB * BYTES_PER_FIELD_ELEMENT;

pub type Bytes32 = [u8; 32];
pub type Bytes48 = [u8; 48];
pub type KZGCommitment = [u8; 48];
pub type KZGProof = [u8; 48];
pub type Blob = [u8; BYTES_PER_BLOB];

#[allow(non_camel_case_types)]
pub type limb_t = u64;

#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_fr {
    l: [limb_t; 256 / 8 / core::mem::size_of::<limb_t>()],
}

#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_fp {
    l: [limb_t; 384 / 8 / core::mem::size_of::<limb_t>()],
}

#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p1 {
    pub x: blst_fp,
    pub y: blst_fp,
    pub z: blst_fp,
}

#[allow(non_camel_case_types)]
#[repr(C)]
pub struct blst_p1_affine {
    pub x: blst_fp,
    pub y: blst_fp,
}

/* 0 is "real" part, 1 is "imaginary" */
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
    todo!()
}

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

    let Ok(poly) = utils::blob_to_polynomial(&input_blob) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(fr_z) = FrElement::from_bytes_be(&z_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let fr_y = poly.evaluate(&fr_z);

    // FIXME: We should not use create_src() for this instantiation.
    let kzg = KZG::new(utils::create_srs());
    let proof = kzg.open(&fr_z, &fr_y, &poly);
    let Ok(compressed_proof): Result<[u8; 48], _> = compress_g1_point(&proof).try_into() else {
        return C_KZG_RET::C_KZG_ERROR; 
    };

    unsafe {
        std::ptr::copy(compressed_proof.as_ptr(), proof_out as *mut u8, BYTES_PER_PROOF);
    }

    C_KZG_RET::C_KZG_OK
}

#[no_mangle]
pub extern "C" fn compute_blob_kzg_proof(
    out: *mut KZGProof,
    blob: *const Blob,
    commitment_bytes: *const Bytes48,
    s: *const KZGSettings,
) -> C_KZG_RET {
    todo!()
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

    let Ok(commitment_g1) = utils::decompress_g1_point(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(z_fr) = FrElement::from_bytes_be(&z_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(y_fr) = FrElement::from_bytes_be(&y_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(proof_g1) = utils::decompress_g1_point(&mut proof_slice) else {
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
    let mut commitment_slice = unsafe { *commitment_bytes };
    let mut proof_slice = unsafe { *proof_bytes };
    let s_struct = unsafe { (*s).clone() };

    let Ok(commitment_g1) = utils::decompress_g1_point(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };
    let Ok(proof_g1) = utils::decompress_g1_point(&mut proof_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    // TODO!!! blob

    todo!()
}

#[no_mangle]
pub extern "C" fn verify_blob_kzg_proof_batch(
    ok: *mut bool,
    blobs: *const Blob,
    commitments_bytes: *const Bytes48,
    proofs_bytes: *const Bytes48,
    n: usize,
    s: *const KZGSettings,
) -> C_KZG_RET {
    todo!()
}
