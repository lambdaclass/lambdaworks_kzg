


pub fn check_point_is_in_subgroup(point: &G1Point) -> bool
{
    let inf = G1Point::neutral_element();
    let aux_point = point.operate_with_self(BLS12381_PRIME_FIELD_ORDER);
    inf == aux_point
}
