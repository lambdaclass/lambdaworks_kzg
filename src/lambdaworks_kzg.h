#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>


/**
 * The number of bytes in a blob.
 */
#define BYTES_PER_BLOB (FIELD_ELEMENTS_PER_BLOB * BYTES_PER_FIELD_ELEMENT)

/**
 * The number of bytes in a KZG commitment.
 */
#define BYTES_PER_COMMITMENT 48

/**
 * The number of bytes in a BLS scalar field element.
 */
#define BYTES_PER_FIELD_ELEMENT 32

/**
 * The number of bytes in a KZG proof.
 */
#define BYTES_PER_PROOF 48

#define FIELD_ELEMENTS_PER_BLOB 4096

/**
 * Order of the subgroup of the curve.
 */
#define TEST_CURVE_1_MAIN_SUBGROUP_ORDER 5

/**
 * Order of the base field (e.g.: order of the coordinates)
 */
#define TEST_CURVE_1_PRIME_FIELD_ORDER 59

typedef enum C_KZG_RET {
  /**
   * Success
   */
  C_KZG_OK = 0,
  /**
   * The supplied data is invalid in some way
   */
  C_KZG_BADARGS,
  /**
   * Internal error - this should never occur
   */
  C_KZG_ERROR,
  /**
   * Could not allocate memory
   */
  C_KZG_MALLOC,
} C_KZG_RET;

typedef uint8_t KZGCommitment[48];

typedef uint8_t Blob[BYTES_PER_BLOB];

typedef uint64_t limb_t;

typedef struct blst_fr {
  limb_t l[4];
} blst_fr;

/**
 *< Internal G2 group element type.
 */
typedef struct blst_fr fr_t;

/**
 * Stores the setup and parameters needed for performing FFTs.
 */
typedef struct FFTSettings {
  /**
   * The maximum size of FFT these settings support, a power of 2.
   */
  uint64_t max_width;
  /**
   *      * We fix a given primitive roots of unity w of order `max_width` via      * `SCALE2_ROOT_OF_UNITY`. Then `expanded_roots_of_unity[i]` == w^i and      * `reverse_roots_of_unity[i]` == w^{-i}. Unusually, both      * `expanded_roots_of_unity` and `reverse_roots_of_unity` have length      * `max_width + 1`. By the above, `expanded_roots_of_unity[max_width] ==      * expanded_roots_of_unity[0] == 1` and similarly for      * `reverse_roots_of_unity`. The redundant element is just there to simplify      * index calculations in some formulas.
   * Ascending powers of the root of unity, length `max_width + 1`.
   */
  fr_t *expanded_roots_of_unity;
  /**
   * Descending powers of the root of unity, length `max_width + 1`.
   */
  fr_t *reverse_roots_of_unity;
  /**
   * Powers of the root of unity in bit-reversal permutation order, length      * `max_width`.
   */
  fr_t *roots_of_unity;
} FFTSettings;

typedef struct blst_fp {
  limb_t l[4];
} blst_fp;

typedef struct blst_p1 {
  struct blst_fp x;
  struct blst_fp y;
  struct blst_fp z;
} blst_p1;

typedef struct blst_p1 g1_t;

typedef struct blst_fp2 {
  struct blst_fp fp[2];
} blst_fp2;

/**
 *< Internal Fr field element type.
 */
typedef struct blst_p2 {
  struct blst_fp2 x;
  struct blst_fp2 y;
  struct blst_fp2 z;
} blst_p2;

/**
 *< Internal G1 group element type.
 */
typedef struct blst_p2 g2_t;

typedef struct KZGSettings {
  /**
   * The corresponding settings for performing FFTs.
   */
  struct FFTSettings *fs;
  /**
   * G1 group elements from the trusted setup,      * in Lagrange form bit-reversal permutation.
   */
  g1_t *g1_values;
  /**
   * G2 group elements from the trusted setup;      * both arrays have `FIELD_ELEMENTS_PER_BLOB` elements.
   */
  g2_t *g2_values;
} KZGSettings;

typedef uint8_t KZGProof[48];

typedef uint8_t Bytes48[48];

typedef uint8_t Bytes32[32];

/**
 * Convert a blob to a KZG commitment.
 *
 * # Params
 *
 * * `out` -  The resulting commitment
 * * `blob` - The blob representing the polynomial to be committed to
 * * `s`    - The trusted setup
 */
enum C_KZG_RET blob_to_kzg_commitment(KZGCommitment *out,
                                      const Blob *blob,
                                      const struct KZGSettings *s);

enum C_KZG_RET compute_blob_kzg_proof(KZGProof *out,
                                      const Blob *blob,
                                      const Bytes48 *commitment_bytes,
                                      const struct KZGSettings *s);

enum C_KZG_RET compute_kzg_proof(KZGProof *proof_out,
                                 Bytes32 *y_out,
                                 const Blob *blob,
                                 const Bytes32 *z_bytes,
                                 const struct KZGSettings *s);

enum C_KZG_RET verify_blob_kzg_proof(bool *ok,
                                     const Blob *blob,
                                     const Bytes48 *commitment_bytes,
                                     const Bytes48 *proof_bytes,
                                     const struct KZGSettings *s);

enum C_KZG_RET verify_blob_kzg_proof_batch(bool *ok,
                                           const Blob *blobs,
                                           const Bytes48 *commitments_bytes,
                                           const Bytes48 *proofs_bytes,
                                           size_t n,
                                           const struct KZGSettings *s);

enum C_KZG_RET verify_kzg_proof(bool *ok,
                                const Bytes48 *commitment_bytes,
                                const Bytes32 *z_bytes,
                                const Bytes32 *y_bytes,
                                const Bytes48 *proof_bytes,
                                const struct KZGSettings *s);
