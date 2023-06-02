use lambdaworks_crypto::commitments::traits::IsCommitmentScheme;
use lambdaworks_kzg::srs::{
    g1_point_to_blst_p1, load_trusted_setup_file,
    load_trusted_setup_file_to_g1_points_and_g2_points, vecs_to_structured_reference_string,
};
use lambdaworks_kzg::traits::Compress;
use lambdaworks_kzg::utils::polynomial_to_blob_with_size;
use lambdaworks_kzg::{
    blst_p1, compute_kzg_proof, free_trusted_setup, kzgsettings_to_structured_reference_string,
    verify_blob_kzg_proof_batch, verify_kzg_proof, Blob, Bytes32, Bytes48, FrElement, G1Point,
    KZGProof, KZGSettings, BYTES_PER_BLOB, C_KZG_RET, FE,
};
use lambdaworks_math::{
    cyclic_group::IsGroup, field::element::FieldElement, polynomial::Polynomial,
    traits::ByteConversion,
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
    let lines = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/trusted_setup.txt"
    ));
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
    let p = G1Point::decompress_g1_point(&mut proof_out).unwrap();
    let inf = G1Point::neutral_element();
    assert_eq!(p, inf);

    let kzg = lambdaworks_kzg::KZG::new(lambdaworks_kzg::utils::create_srs());
    let commitment = kzg.commit(&polynomial);
    let commitment_bytes = G1Point::compress_g1_point(&commitment).unwrap();

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
    let lines = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/trusted_setup.txt"
    ));
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
    let commitment_bytes = G1Point::compress_g1_point(&commitment).unwrap();

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
    let proof_out_point = G1Point::decompress_g1_point(&mut proof_out).unwrap();
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
    let ok_enum_kzg = C_KZG_RET::C_KZG_OK;

    // read srs from file
    let lines = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/trusted_setup.txt"
    ));
    let lines = lines.lines();
    let s = load_trusted_setup_file(lines).unwrap();

    // assign poly 1
    let polynomial = Polynomial::<FrElement>::new(&[FieldElement::one()]);
    let blob = polynomial_to_blob_with_size(&polynomial).unwrap();
    // z = 1
    let z_bytes: Bytes32 = FE::from(1).to_bytes_be().try_into().unwrap();

    // assign poly2 = x
    let polynomial2 = Polynomial::<FrElement>::new(&[FieldElement::zero(), FieldElement::one()]);
    let blob2 = polynomial_to_blob_with_size(&polynomial2).unwrap();
    // z = 2
    let z_bytes2: Bytes32 = FE::from(2).to_bytes_be().try_into().unwrap();
    let mut blob_both: [u8; BYTES_PER_BLOB * 2] = [0; BYTES_PER_BLOB * 2];
    blob_both[..BYTES_PER_BLOB].copy_from_slice(&blob);
    blob_both[BYTES_PER_BLOB..].copy_from_slice(&blob2);

    // output buffers 1
    let mut proof_out1: KZGProof = [0; 48];
    let mut y_out1: Bytes32 = [0; 32];

    // output buffers 2
    let mut proof_out2: KZGProof = [0; 48];
    let mut y_out2: Bytes32 = [0; 32];

    let ret = compute_kzg_proof(
        &mut proof_out1 as *mut KZGProof,
        &mut y_out1 as *mut Bytes32,
        &blob as *const Blob,
        &z_bytes as *const Bytes32,
        &s as *const KZGSettings,
    );

    // assert ret function
    assert_eq!(ret, ok_enum_kzg);

    let ret = compute_kzg_proof(
        &mut proof_out2 as *mut KZGProof,
        &mut y_out2 as *mut Bytes32,
        &blob2 as *const Blob,
        &z_bytes2 as *const Bytes32,
        &s as *const KZGSettings,
    );

    // assert ret function
    assert_eq!(ret, ok_enum_kzg);

    // build proof_batch
    let mut proof_out_batch = [0_u8; 96];
    proof_out_batch[..48].copy_from_slice(&proof_out1);
    proof_out_batch[48..].copy_from_slice(&proof_out2);

    // build y_out_batch
    let mut y_out_batch = [0_u8; 64];
    y_out_batch[..32].copy_from_slice(&y_out1);
    y_out_batch[32..].copy_from_slice(&y_out2);

    let srs = kzgsettings_to_structured_reference_string(&s).unwrap();
    let kzg = lambdaworks_kzg::KZG::new(srs);
    let commitment1 = kzg.commit(&polynomial);
    let commitment_bytes1 = G1Point::compress_g1_point(&commitment1).unwrap();

    let commitment2 = kzg.commit(&polynomial2);
    let commitment_bytes2 = G1Point::compress_g1_point(&commitment2).unwrap();
    let mut commitment_bytes_batch = [0_u8; 96];
    commitment_bytes_batch[..48].copy_from_slice(&commitment_bytes1);
    commitment_bytes_batch[48..].copy_from_slice(&commitment_bytes2);

    // verify blob as a batch
    let mut ok = false;
    verify_blob_kzg_proof_batch(
        &mut ok as *mut bool,
        blob_both.as_ptr() as *const Blob,
        commitment_bytes_batch.as_ptr() as *const Bytes48,
        proof_out_batch.as_ptr() as *const Bytes48,
        2,
        &s as *const KZGSettings,
    );

    assert!(ok);
    assert_eq!(ret, ok_enum_kzg);
}

#[test]
fn test_read_srs() {
    let lines = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/trusted_setup.txt"
    ));
    let lines = lines.lines();
    let (g1_points, g2_points) = load_trusted_setup_file_to_g1_points_and_g2_points(lines).unwrap();

    let g = g1_points[0].clone();
    let compressed = G1Point::compress_g1_point(&g).unwrap();
    let hex_string = hex::encode(compressed);

    assert_eq!(hex_string,
            "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");

    let srs_from_file = vecs_to_structured_reference_string(&g1_points, &g2_points);
    let lines = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/trusted_setup.txt"
    ));
    let lines = lines.lines();
    let mut s = load_trusted_setup_file(lines).unwrap();
    let srs = kzgsettings_to_structured_reference_string(&s).unwrap();

    assert_eq!(srs_from_file, srs);

    // free memory used by srs struct
    unsafe { free_trusted_setup(&mut s as *mut KZGSettings) };
}
