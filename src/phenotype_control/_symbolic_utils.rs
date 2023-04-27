use biodivine_lib_bdd::{Bdd, BddPartialValuation, BddVariable, BddVariableSet};

/// Build a BDD which is true for all valuations of the given `variables` that contain up to
/// `bound` true values.
///
/// For example, given `bound=0`, the result is a single conjunctive clause stating that all
/// variables must be false. If `bound=1`, the resulting BDD is satisfied for any valuation
/// of `variables` where at most one variable is satisfied, and so on...
pub fn mk_bdd_up_to_bound(ctx: &BddVariableSet, variables: &[BddVariable], bound: usize) -> Bdd {
    let mut partial_valuation = BddPartialValuation::empty();
    for var in variables {
        partial_valuation.set_value(*var, false);
    }
    let mut result = ctx.mk_conjunctive_clause(&partial_valuation);
    for _i in 0..bound {
        let mut result_plus_one = result.clone();
        for var in variables {
            // First, filter valuations in `result` to those which do not include
            // `var`. Then flip the value such that the new BDD actually includes valuations
            // with `var`. In other words, take valuations of size `k` without `var`,
            // and then add `var` to all of them, resulting in valuations of size `k+1`.
            let not_var = ctx.mk_literal(*var, false);
            let with_var = Bdd::fused_binary_flip_op(
                (&result, None),
                (&not_var, None),
                Some(*var),
                biodivine_lib_bdd::op_function::and
            );
            result_plus_one = result_plus_one.or(&with_var);
        }
        result = result_plus_one;
    }

    result
}

/// The same as `mk_bdd_up_to_bound`, but the result only allows
/// valuations of size exactly `bound`.
pub fn mk_bdd_of_bound(ctx: &BddVariableSet, variables: &[BddVariable], bound: usize) -> Bdd {
    let mut partial_valuation = BddPartialValuation::empty();
    for var in variables {
        partial_valuation.set_value(*var, false);
    }
    let mut result = ctx.mk_conjunctive_clause(&partial_valuation);
    for _i in 0..bound {
        let mut result_plus_one = ctx.mk_false();
        for var in variables {
            // First, filter valuations in `result` to those which do not include
            // `var`. Then flip the value such that the new BDD actually includes valuations
            // with `var`. In other words, take valuations of size `k` without `var`,
            // and then add `var` to all of them, resulting in valuations of size `k+1`.
            let not_var = ctx.mk_literal(*var, false);
            let with_var = Bdd::fused_binary_flip_op(
                (&result, None),
                (&not_var, None),
                Some(*var),
                biodivine_lib_bdd::op_function::and
            );
            result_plus_one = result_plus_one.or(&with_var);
        }
        result = result_plus_one;
    }

    result
}

#[cfg(test)]
mod tests {
    use biodivine_lib_bdd::BddVariableSet;
    use crate::phenotype_control::_symbolic_utils::{mk_bdd_of_bound, mk_bdd_up_to_bound};

    #[test]
    pub fn test_bounded_bdd() {
        let vars = BddVariableSet::new_anonymous(5);

        let bdd = mk_bdd_of_bound(&vars, &vars.variables(), 3);
        // The number of such valuations is exactly the binomial coefficient.
        assert_eq!(bdd.cardinality(), binomial(5, 3) as f64);
    }

    #[test]
    pub fn test_bound_up_to_bdd() {
        // The same as above, but the result is a sum of binomials.
        let vars = BddVariableSet::new_anonymous(5);

        let bdd = mk_bdd_up_to_bound(&vars, &vars.variables(), 3);
        // The number of such valuations is exactly the binomial coefficient.
        assert_eq!(bdd.cardinality(), (binomial(5, 3) + binomial(5, 2) + binomial(5, 1) + binomial(5, 0)) as f64);
    }

    pub fn binomial(n: usize, k: usize) -> usize {
        factorial(n) / (factorial(k) * factorial(n - k))
    }

    pub fn factorial(x: usize) -> usize {
        if x == 0 {
            1
        } else {
            x * factorial(x - 1)
        }
    }

}