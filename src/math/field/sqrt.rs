use crate::math::unsigned_integer::traits::IsUnsignedInteger;

use super::traits::IsField;

/// Represents whether an integer is a quadratic residue, non-residue (non-zero) or zero.
/// Discriminants are based on Euler's criterion.
pub enum Quadratic {
    NonResidue = -1,
    Zero = 0,
    NonZeroResidue = 1,
}

/// Gives a field the ability of calculating an element's legendre symbol and to determine
/// quadratic residuity.
pub trait LegendreSymbol<F>
where
    F: IsField,
    F::BaseType: IsUnsignedInteger,
{
    const MODULUS_MINUS_ONE_DIV_TWO: F::BaseType;

    /// Returns `a` Legendre's symbol.
    fn legendre_symbol(a: &F::BaseType) -> F::BaseType {
        F::pow(a, Self::MODULUS_MINUS_ONE_DIV_TWO)
    }

    /// Returns a `Quadratic` enum representing `a`'s quadratic residuity.
    fn euler_criterion(a: &F::BaseType) -> Quadratic {
        match Self::legendre_symbol(a) {
            x if x == F::neg(&F::one()) => Quadratic::NonResidue,
            x if x == F::zero() => Quadratic::Zero,
            x if x == F::one() => Quadratic::NonZeroResidue,
            _ => unreachable!(),
        }
    }
}
