use crate::control::ControlMap;
use biodivine_lib_bdd::Bdd;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColoredVertices;
use biodivine_lib_param_bn::VariableId;

impl ControlMap {
    /// Remove from this control map any results that *do not* perturb `variable`.
    /// If `value` is given, only keep perturbations which result in this value.
    pub fn require_perturbation(&mut self, variable: VariableId, value: Option<bool>) {
        let require = self.context.fix_perturbation(variable, value);
        self.perturbation_set = self.perturbation_set.intersect(&require);
    }

    /// Remove from this control map any results that perturb `variable`. If `value` is given,
    /// only remove perturbations which result in this value.
    pub fn exclude_perturbation(&mut self, variable: VariableId, value: Option<bool>) {
        let exclude = self.context.fix_perturbation(variable, value);
        self.perturbation_set = self.perturbation_set.minus(&exclude);
    }

    pub fn as_bdd(&self) -> &Bdd {
        self.perturbation_set.as_bdd()
    }

    pub fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

    pub fn controllable_colors(&self) -> Bdd {
        let bdd_context = self.context.as_symbolic_context();
        let mut bdd = self.perturbation_set.colors().into_bdd();
        for v in self.context.variables() {
            let parameter = self.context.get_perturbation_parameter(v);
            if let Some(parameter) = parameter {
                for (_, bdd_var) in bdd_context.get_explicit_function_table(parameter) {
                    // There should be only one because arity is zero
                    bdd = bdd.var_exists(bdd_var);
                }
            }
        }

        // Now just fix the control parameters and variables to true so that the final cardinality is correct.
        for s_var in bdd_context.state_variables() {
            bdd = bdd.var_select(*s_var, true);
        }
        for v in self.context.variables() {
            let parameter = self.context.get_perturbation_parameter(v);
            if let Some(parameter) = parameter {
                for (_, bdd_var) in bdd_context.get_explicit_function_table(parameter) {
                    // There should be only one because arity is zero
                    bdd = bdd.var_select(bdd_var, true);
                }
            }
        }

        bdd
    }

    /// Compute the set of original colors (without perturbations) that can be controlled
    /// using this control map.
    pub fn controllable_colors_cardinality(&self) -> f64 {
        self.controllable_colors().cardinality()
    }

    /// Compute the number of vertices the source can jump to due to different perturbations.
    pub fn jump_vertices(&self) -> f64 {
        self.perturbation_set.vertices().approx_cardinality()
    }
}
