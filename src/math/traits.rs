use crate::math::{
    elliptic_curve::short_weierstrass::errors::DeserializationError,
    errors::ByteConversionError,
    field::{element::FieldElement, traits::IsField},
};

/// A trait for converting an element to and from its byte representation and
/// for getting an element from its byte representation in big-endian or
/// little-endian order.
pub trait ByteConversion {
    /// Returns the byte representation of the element in big-endian order.
    fn to_bytes_be(&self) -> Vec<u8>;

    /// Returns the byte representation of the element in little-endian order.
    fn to_bytes_le(&self) -> Vec<u8>;

    /// Returns the element from its byte representation in big-endian order.
    fn from_bytes_be(bytes: &[u8]) -> Result<Self, ByteConversionError>
    where
        Self: std::marker::Sized;

    /// Returns the element from its byte representation in little-endian order.
    fn from_bytes_le(bytes: &[u8]) -> Result<Self, ByteConversionError>
    where
        Self: std::marker::Sized;
}

/// Serialize function without args
/// Used for serialization when formatting options are not relevant
pub trait Serializable {
    /// Default serialize without args
    fn serialize(&self) -> Vec<u8>;
}

/// Deserialize function without args
/// Used along with the Serializable trait
pub trait Deserializable {
    fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError>
    where
        Self: Sized;
}

pub trait IsRandomFieldElementGenerator<F: IsField> {
    fn generate(&self) -> FieldElement<F>;
}
