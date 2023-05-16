use criterion::*;
use lambdaworks_kzg::{
    srs::kzgsettings_to_structured_reference_string, srs::load_trusted_setup_file,
};

fn bench_read_srs(c: &mut Criterion) {
    const TRUSTED_SETUP_FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/test/trusted_setup.txt");

    let s = load_trusted_setup_file(TRUSTED_SETUP_FILE.lines()).unwrap();
    c.bench_function("kzgsettings_to_ref_string", |b| {
        b.iter(|| black_box(kzgsettings_to_structured_reference_string(black_box(&s))));
    });
}

criterion_group!(lambdaworks, bench_read_srs);
criterion_main!(lambdaworks);
