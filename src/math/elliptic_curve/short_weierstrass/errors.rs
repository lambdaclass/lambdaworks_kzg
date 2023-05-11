use thiserror::Error;

use crate::math::errors::ByteConversionError;

#[derive(Error, Debug, PartialEq, Eq)]
pub enum DeserializationError {
    #[error("Invalid amount of bytes")]
    InvalidAmountOfBytes,
    #[error("Error when creating a field from bytes")]
    FieldFromBytesError,
    #[error("Error trying to load a pointer bigger than the supported architecture")]
    PointerSizeError,
}

impl From<ByteConversionError> for DeserializationError {
    fn from(error: ByteConversionError) -> Self {
        match error {
            ByteConversionError::FromBEBytesError => DeserializationError::FieldFromBytesError,
            ByteConversionError::FromLEBytesError => DeserializationError::FieldFromBytesError,
            ByteConversionError::InvalidValue => todo!(),
            ByteConversionError::PointNotInSubgroup => todo!(),
            ByteConversionError::ValueNotCompressed => todo!(),
        }
    }
}
