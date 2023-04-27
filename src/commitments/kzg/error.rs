use std::fmt::Display;

#[derive(Debug)]
pub enum KzgError {
    BadArgs,
}

impl Display for KzgError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            KzgError::BadArgs => write!(f, "The supplied data is invalid"),
        }
    }
}
