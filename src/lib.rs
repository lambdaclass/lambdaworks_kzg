#![allow(unused_variables)]

pub mod commitments;
pub mod math;
pub mod utils;

use commitments::kzg::{FrElement, StructuredReferenceString, G1};
use math::{
    elliptic_curve::{
        short_weierstrass::{
            curves::bls12_381::{
                curve::BLS12381Curve, field_extension::BLS12381PrimeField,
                pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
            },
            point::ShortWeierstrassProjectivePoint,
        },
        traits::IsPairing,
    },
    errors::ByteConversionError,
    field::element::FieldElement,
    polynomial::Polynomial,
    traits::ByteConversion,
    unsigned_integer::element::U256,
};
use rand::Rng;
use std::{marker, mem::size_of, slice};

use crate::math::{
    cyclic_group::IsGroup,
    elliptic_curve::traits::{FromAffine, IsEllipticCurve},
};

pub type G1Point = ShortWeierstrassProjectivePoint<BLS12381Curve>;
pub type G2Point = ShortWeierstrassProjectivePoint<BLS12381TwistCurve>;

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
    todo!()
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

    let Ok(commitment_g1) = get_point_from_bytes(&mut commitment_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(z_fr) = FrElement::from_bytes_be(&z_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(y_fr) = FrElement::from_bytes_be(&y_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    let Ok(proof_g1) = get_point_from_bytes(&mut proof_slice) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

    // verify pairing
    let Ok(ret) = verify_kzg_proof_impl(&commitment_g1, &z_fr, &y_fr, &proof_g1, &s_struct) else {
        return C_KZG_RET::C_KZG_ERROR;
    };

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

#[allow(dead_code)]
fn blob_to_polynomial(
    blob: *const Blob,
) -> Result<Polynomial<FrElement>, crate::math::errors::ByteConversionError>
where
    FrElement: ByteConversion,
{
    let input_blob = unsafe { *blob };
    let mut coefficients = Vec::new();

    for elem_bytes in input_blob.chunks(BYTES_PER_FIELD_ELEMENT) {
        let f = FrElement::from_bytes_le(elem_bytes)?;
        coefficients.push(f);
    }

    Ok(Polynomial::new(&coefficients))
}

fn get_point_from_bytes(input_bytes: &mut [u8; 48]) -> Result<G1Point, ByteConversionError> {
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
    // TODO: Use optimized sqrt function
    let y = y_squared.inv().pow(2_u16);

    let point = G1Point::from_affine(x, y).map_err(|_| ByteConversionError::InvalidValue)?;

    if utils::check_point_is_in_subgroup(&point) {
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
    //todo!()
}

fn verify_kzg_proof_impl(
    commitment: &G1Point,
    z_fr: &FrElement,
    y_fr: &FrElement,
    proof_g1: &G1Point,
    s_struct: &KZGSettings,
) -> Result<bool, C_KZG_RET> {
    let src = create_srs();

    // let g2_gen = BLS12381TwistCurve::generator();
    // let x_g2 = g2_gen.operate_with_self(z_fr.representative());

    // let s_g2_values =
    //     unsafe { slice::from_raw_parts(s_struct.g2_values, std::mem::size_of::<g2_t>()) };

    // let s_g2_values_1 = s_g2_values[1];

    // let X_minus_z = s_g2_values - x_g2;

    todo!()
}

fn create_srs() -> StructuredReferenceString<
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
mod test {
    use super::*;
}
