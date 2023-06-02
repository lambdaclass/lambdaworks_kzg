use lambdaworks_math::{cyclic_group::IsGroup, elliptic_curve::traits::EllipticCurveError};

pub trait Compress {
    type G1Point: IsGroup;
    type G2Point: IsGroup;

    fn compress_g1_point(point: &Self::G1Point) -> Result<[u8; 48], EllipticCurveError>;

    fn decompress_g1_point(input_bytes: &mut [u8; 48])
        -> Result<Self::G1Point, EllipticCurveError>;

    fn decompress_g2_point(input_bytes: &mut [u8; 96])
        -> Result<Self::G2Point, EllipticCurveError>;
}
