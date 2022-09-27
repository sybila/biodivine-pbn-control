from biodivine_aeon import *
from pathlib import Path
import sys
import random
import math 
import os

# Get a list of variable names appearing in the given network.
def network_variable_names(network):
	return list(map(lambda id: network.get_variable_name(id), network.variables()))

# Make a copy of the given BN where all regulations have observability removed.
def clear_observability(network):
	variables = network_variable_names(network)
	new_rg = RegulatoryGraph(variables)
	
	for regulation in network.graph().regulations():
		regulation['observable'] = False
		regulation['source'] = network.get_variable_name(regulation['source'])
		regulation['target'] = network.get_variable_name(regulation['target'])
		new_rg.add_regulation(regulation)

	new_bn = BooleanNetwork(new_rg)

	for variable in variables:		
		new_bn.set_update_function(variable, network.get_update_function(variable))

	return new_bn

# Produce a new network where all the given variables have erased update functions.
def erase_functions(network, erase):
	result = BooleanNetwork(network.graph())
	for variable in network.variables():
		variable = network.get_variable_name(variable)
		if variable not in erase:
			result.set_update_function(variable, network.get_update_function(variable))		
	return result


# Group network varaibles based on their arity.
def get_arity_list(network):
	rg = network.graph()
	variables = network_variable_names(network)

	def push_arity(list, name, arity):
		while len(list) <= arity:
			list.append([])
		list[arity].append(name)

	arities = []
	for variable in variables:
		arity = len(rg.regulators(variable))
		push_arity(arities, variable, arity)

	return arities

# Compute a list of states that appear within the attractors of the given network.
# For coloured networks, not all attractors are necessarily covered.
def sample_attractors(network):
	stg = SymbolicAsyncGraph(network)
	attractors = find_attractors(stg)

	sampled_states = []
	for attractor_set in attractors:
		vertex = attractor_set.pick_vertex().vertices().list_vertices()[0]
		sampled_states.append(vertex)

	return sampled_states

# Return the number of colors in the given network for which the seed state
# is an attractor. If the computation exceeds the given complexity bound,
# it is terminated and `None` is returned.
def count_attractor_colors(network, seed, bound=1_000_000):
	stg = SymbolicAsyncGraph(network)
	seed = stg.fix_vertex(seed)
	
	bwd = seed
	while True:
		done = True
		for var in reversed(model.variables()):
			step = stg.var_pre(var, bwd)
			if not step.is_subset(bwd):
				bwd = bwd.union(step)
				done = False
				break
		if done:
			break
		if bwd.symbolic_size() > bound:
			return None
	
	attractor = seed
	while True:
		done = True
		for var in reversed(network.variables()):		
			step = stg.var_post(var, attractor)

			if not step.is_subset(attractor):
				bad_colors = step.minus(bwd).colors()
				attractor = attractor.union(step).minus_colors(bad_colors)
				done = False
				break
		if done:
			return attractor.colors().cardinality()
		if attractor.symbolic_size() > bound:
			return None

# Return a random value from the given list of seeds that is different
# than the given seed.
def other_random_seed(seeds, seed):
	assert len(seeds) > 1
	while True:		
		# I know this is inefficient and I don't care.
		chosen = random.choice(seeds)				
		if chosen != seed:
			return chosen





# <<<<<<<<<<<<<<<<<<<<<< MAIN >>>>>>>>>>>>>>>>>>>>>>>>>





BASE_MODEL_PATH = sys.argv[1]
MAX_ERASED_ARITY = 5
MAX_ERASED_COUNT = 10
SAMPLES_PER_ERASE_COUNT = 20
EXPERIMENTS_PER_GROUP = 20
random.seed(1234567890)

# Prepare directory for writing outputs
BASE_MODEL_FILE_NAME = os.path.basename(BASE_MODEL_PATH)
BASE_MODEL_NAME = BASE_MODEL_FILE_NAME.replace(".aeon", "").replace("_witness", "")
os.mkdir(BASE_MODEL_NAME)

model = BooleanNetwork.from_aeon(Path(BASE_MODEL_PATH).read_text())

# PerturbedGraph does not respect observability constraints, hence we remove them
# in order to get comparable results.
model = clear_observability(model)

attractor_seeds = sample_attractors(model)

model_rg = model.graph()
erasure_candidates = list(filter(lambda var: len(model_rg.regulators(var)) <= MAX_ERASED_ARITY, network_variable_names(model)))

print("Erasure candidates:", len(erasure_candidates))

groups = []

def append_group(groups, magnitude, item):
	while len(groups) <= magnitude:
		groups.append([])
	groups[magnitude].append(item)	

for fn_count in range(1,min(MAX_ERASED_COUNT, len(erasure_candidates))+1):
	print(f">>>>>>>>>>>>>>>>>>>>> ERASE {fn_count} <<<<<<<<<<<<<<<<<<<<<<<<<")
	for i in range(0, SAMPLES_PER_ERASE_COUNT):
		random.shuffle(erasure_candidates)

		erased = erasure_candidates[0:fn_count]
		erased_network = erase_functions(model, erased)
		print("Total colors:", round(math.log2(SymbolicAsyncGraph(erased_network).unit_colors().cardinality())))

		for seed in attractor_seeds:
			# Check if this combination of seed and erased variables already exists somewhere.
			skip = False
			for group in groups:
				for existing in group:
					if existing[0] == seed and set(existing[1]) == set(erased):
						skip = True
			# If it does, just skip it.						
			if skip:
				print("-", end=" ", flush=True)
				continue

			count = count_attractor_colors(erased_network, seed)
			if count != None and count > 0:
				magnitude = math.floor(math.log2(count) / 5)
				append_group(groups, magnitude, (seed, erased, erased_network))
				print(f"{round(math.log2(count))}/{magnitude}", end=" ", flush=True)				
			else:
				print("x", end=" ", flush=True)
		print("")


for i in range(0, len(groups)):
	print(f"Group [2^{i*5} ... 2^{(i+1)*5}]: {len(groups[i])}")
	out_dir = f"{BASE_MODEL_NAME}/group_{i+1}"	
	os.mkdir(out_dir)
	group = groups[i]
	random.shuffle(group)
	for j in range(0, min(EXPERIMENTS_PER_GROUP, len(group))):		
		out_file = f"{out_dir}/{j+1}.aeon"		
		(target, erased, erased_network) = group[j]
		source = other_random_seed(attractor_seeds, target)
		out_content = f"#source:{source}\n#target:{target}\n{erased_network.to_aeon()}"
		Path(out_file).write_text(out_content)











