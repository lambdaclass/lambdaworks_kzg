use lambdaworks_kzg::commitments::traits::IsCommitmentScheme;
use lambdaworks_kzg::compress::{compress_g1_point, decompress_g1_point};
use lambdaworks_kzg::math::{
    cyclic_group::IsGroup, field::element::FieldElement, polynomial::Polynomial,
    traits::ByteConversion,
};
use lambdaworks_kzg::srs::{
    g1_point_to_blst_p1, load_trusted_setup_file,
    load_trusted_setup_file_to_g1_points_and_g2_points, vecs_to_structured_reference_string,
};
use lambdaworks_kzg::utils::polynomial_to_blob_with_size;
use lambdaworks_kzg::{
    blst_p1, compute_kzg_proof, free_trusted_setup, kzgsettings_to_structured_reference_string,
    verify_kzg_proof, Blob, Bytes32, Bytes48, FrElement, G1Point, KZGProof, KZGSettings, C_KZG_RET,
    FE,
};
use pretty_assertions::assert_eq;

#[test]
fn test_compute_kzg_proof_for_a_simple_poly() {
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

    // read srs from file
    let lines = std::fs::read_to_string("tests/trusted_setup.txt").unwrap();
    let lines = lines.lines();
    let mut s = load_trusted_setup_file(lines).unwrap();

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

    let kzg = lambdaworks_kzg::KZG::new(lambdaworks_kzg::utils::create_srs());
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

    assert!(ok);
    assert_eq!(ret_verify, ok_enum_kzg);

    // free memory used by srs struct
    unsafe { free_trusted_setup(&mut s as *mut KZGSettings) };
}

#[test]
fn test_compute_kzg_proof_for_a_simple_poly_2() {
    // Test this case:
    // polinomial: x
    // evaluation fr: 2
    // Expected result:
    // proof: first point of srs
    // y_out: 2

    // output buffers
    let mut proof_out: KZGProof = [0; 48];
    let mut y_out: Bytes32 = [0; 32];

    // assign poly = x
    let polynomial = Polynomial::<FrElement>::new(&[FieldElement::zero(), FieldElement::one()]);
    let blob = polynomial_to_blob_with_size(&polynomial).unwrap();
    // z = 2
    let z_bytes: Bytes32 = FE::from(2).to_bytes_be().try_into().unwrap();
    let z_fr = FrElement::from_bytes_be(&z_bytes).unwrap();

    // read srs from file
    let lines = std::fs::read_to_string("tests/trusted_setup.txt").unwrap();
    let lines = lines.lines();
    let mut s = load_trusted_setup_file(lines).unwrap();

    let ret = compute_kzg_proof(
        &mut proof_out as *mut KZGProof,
        &mut y_out as *mut Bytes32,
        &blob as *const Blob,
        &z_bytes as *const Bytes32,
        &s as *const KZGSettings,
    );

    let y_out_field = FrElement::from_bytes_be(&y_out).unwrap();
    assert_eq!(y_out_field, z_fr);

    // assert ret function
    let ok_enum_kzg = C_KZG_RET::C_KZG_OK;
    assert_eq!(ret, ok_enum_kzg);

    // verify proof
    let srs = kzgsettings_to_structured_reference_string(&s).unwrap();
    let kzg = lambdaworks_kzg::KZG::new(srs);
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

    assert!(ok);
    assert_eq!(ret_verify, ok_enum_kzg);

    // get first point of the srs
    // assert proof is the first element in g2 points of srs
    let g1_values_slice: &[blst_p1] =
        unsafe { std::slice::from_raw_parts(s.g1_values as *const blst_p1, 4096) };
    let first_point_srs_blst = g1_values_slice[0];
    let proof_out_point = decompress_g1_point(&mut proof_out).unwrap();
    let proof_out_point_blst = g1_point_to_blst_p1(&proof_out_point);
    assert_eq!(first_point_srs_blst, proof_out_point_blst);

    // check commitment is second point of SRS
    let second_point_srs_blst = g1_values_slice[1];
    let commitment_blst = g1_point_to_blst_p1(&commitment);
    assert_eq!(second_point_srs_blst, commitment_blst);

    // free memory used by srs struct
    unsafe { free_trusted_setup(&mut s as *mut KZGSettings) };
}

#[test]
fn test_batch_proof() {
    /*
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
    */
}

#[test]
fn test_read_srs() {
    let lines = std::fs::read_to_string("tests/trusted_setup.txt").unwrap();
    let lines = lines.lines();
    let (g1_points, g2_points) = load_trusted_setup_file_to_g1_points_and_g2_points(lines).unwrap();

    let g = g1_points[0].clone();
    let compressed = compress_g1_point(&g).unwrap();
    let hex_string = hex::encode(compressed);

    assert_eq!(hex_string,
            "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");

    let srs_from_file = vecs_to_structured_reference_string(&g1_points, &g2_points);
    let lines = std::fs::read_to_string("tests/trusted_setup.txt").unwrap();
    let lines = lines.lines();
    let mut s = load_trusted_setup_file(lines).unwrap();
    let srs = kzgsettings_to_structured_reference_string(&s).unwrap();

    assert_eq!(srs_from_file, srs);

    // free memory used by srs struct
    unsafe { free_trusted_setup(&mut s as *mut KZGSettings) };
}
