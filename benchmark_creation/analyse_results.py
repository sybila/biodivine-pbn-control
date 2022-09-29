# This script is for reading the outputs of experiments and extracting useful statistics from them.

from pathlib import Path
import os
import sys
import re

# Regexes for parsing metadata from benchmark output.
RE_PERTURBATION_COLORS = re.compile("Perturbation colors: (\\d+)")
RE_ATTRACTOR_COLORS = re.compile("Attractor colors: (\\d+)")
RE_ELAPSED_MS = re.compile("Elapsed: (\\d+) ms")
RE_MIN_ARITY = re.compile("Lowest cardinality: (\\d+)")
RE_MAX_ARITY = re.compile("Highest cardinality: (\\d+)")
RE_NUM_FUNCTIONS = re.compile("Unknown update functions: (\\d+)")


RESULTS_DIR = sys.argv[1]
#print(f"Reading benchmarks from {RESULTS_DIR}")

print("Benchmark name, Attractor colors, Unknown functions, Arity range, Elapsed (ms)")

for result in os.listdir(RESULTS_DIR):
	if not result.endswith("_out.txt"):
		continue

	result_text = Path(f"{RESULTS_DIR}/{result}").read_text()

	name = result.replace("_out.txt", "")
	if RE_ELAPSED_MS.search(result_text):		
		elapsed = RE_ELAPSED_MS.search(result_text).group(1)
		perturbation_colors = RE_PERTURBATION_COLORS.search(result_text).group(1)
		attractor_colors = RE_ATTRACTOR_COLORS.search(result_text).group(1)
		min_arity = RE_MIN_ARITY.search(result_text).group(1)
		max_arity = RE_MAX_ARITY.search(result_text).group(1)
		num_functions = RE_NUM_FUNCTIONS.search(result_text).group(1)


		colors = int(int(attractor_colors) / int(perturbation_colors))
		print(f"{name}, {colors}, {num_functions}, {min_arity}-{max_arity}, {elapsed}")
	else:
		print(f"{name}, x, x, x, x")