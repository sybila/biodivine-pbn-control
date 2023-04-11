# This is a modified version of https://github.com/daemontus/artefact-aeon-py/blob/2f59d32acc8ac278f18a4cf6fa30e6b4b3c0c619/run.py
# In particular, the changes are:
#  - We don't invoke UNIX `time` explicitly. Instead we assume that the binary will print the time in a format that we can understand.
#  - The benchmarks are taken from a two-level hierarchy of folders instead of a single directory (this influences the output as well).
#  -
import sys
import os
import re
import time
from multiprocessing import Process

# Spawn a new external process.
def SPAWN(command):
    print("Spawn:", command)
    result = os.system(command)
    sys.exit(result)

# Regex used for parsing the standardised output of
# UNIX time utility (when invoked with -p argument)
RE_TIME = re.compile("\\s*real\\s*(\\d+\\.?\\d*)\\s*")

models = [
    ("cardiac", "FHF"),
    ("cardiac", "SHF"),
    ("cardiac", "None"),
    ("reduced_mapk", "apoptosis"),
    ("reduced_mapk", "proliferation"),
    ("reduced_mapk", "no_decision"),
    ("reduced_mapk", "growth_arrest"),
    ("erbb", "phosporilation"),
    ("erbb", "non_phosporilation"),
    ("tumour", "apoptosis"),
    ("tumour", "emt"),
    ("tumour", "hybrid"),
    ("tumour", "metastasis"),
    ("cell_fate", "apoptosis"),
    ("cell_fate", "survival"),
    ("cell_fate", "naive"),
    ("cell_fate", "necrosis"),
    ("full_mapk", "apoptosis"),
    ("full_mapk", "proliferation"),
    ("full_mapk", "no_decision"),
    ("full_mapk", "growth_arrest"),
    ("t_lgl", "apoptosis"),
    ("t_lgl", "proliferation")
]

if __name__ == "__main__":
    print(">>>>>>>>>> START BENCHMARK RUN")

    # 2 days
    CUT_OFF = "48h"
    SCRIPT = "./target/release/experiment_scalability"
    INTERACTIVE = False
    PARALLEL = 4

    PERTURBATION_MAX_SIZE = "3"

    # Set timeout binary based on OS (macOS needs gtimeout)
    TIMEOUT = 'none'

    if TIMEOUT == 'none':
        code = os.system('timeout --help > /dev/null 2>&1')
        if code == 0:
            TIMEOUT = 'timeout'
            print("Timeout utility ok.")

    if TIMEOUT == 'none':
        code = os.system('gtimeout --help > /dev/null 2>&1')
        if code == 0:
            TIMEOUT = 'gtimeout'
            print("Timeout utility ok.")

    if TIMEOUT == 'none':
        print('ERROR: No timeout utility found.')
        exit()


    # Create output directory
    OUT_DIR = "phenotype_output_bench"
    OUT_DIR = "_run_" + OUT_DIR + "_" + str(int(time.time()))
    os.mkdir(OUT_DIR)

    # Create output stats file
    TIMES = open(OUT_DIR + "/"+ "phenotypes" + "_times.csv", "w")
    TIMES.write("Benchmark, Time[s]\n")

    # Create an aggregated stats file
    AGGREGATION = open(OUT_DIR + "/"+ "phenotypes" + "_aggregated.csv", "w")
    AGGREGATION.write("Time[s], No. Completed\n")

    # Here, save all runtimes.
    AGGREGATION_LIST = []

    BENCHMARKS = [i for i in reversed(range(1,52))]

    MAX_MEM = {}
    # Handle data from a finished process. In particular,
    # update AGGREGATION_LIST and TIMES file.
    def PROCESS_RESULT(process, name, output_file):
        print("Finished:", output_file)
        # print(pid)
        # print("Max mem: ", MAX_MEM.get(pid, -1))
        is_success = process.exitcode == 0
        with open(output_file, 'r') as f:
            lines = f.read().splitlines()
            # Try to parse runtime statistics.
            # Time stats are three lines from the end.
            if len(lines) >= 3:
                time_line = lines[-3]
                if RE_TIME.match(time_line) and is_success:
                    # Success, we found time!
                    time = str(RE_TIME.match(time_line).group(1))
                    print("Success. Elapsed: ", time)
                    TIMES.write(name + ", " + time + "\n")
                    AGGREGATION_LIST.append(float(time))
                else:
                    # Fail: output exists but does not have
                    # correct time format.
                    print("Fail. Last line of output:")
                    print(lines[-1])
                    TIMES.write(name + ", " + "fail" + "\n")
            elif len(lines) > 0:
                # Fail: There is some output, but not enough
                # for a successful process.
                print("Fail. Last line of output:")
                print(lines[-1])
                TIMES.write(name + ", " + "fail" + "\n")
            else:
                # Fail: There is no output.
                print("Fail. No output found.")
                TIMES.write(name + ", " + "fail" + "\n")
        TIMES.flush()

    ACTIVE = []
    while len(ACTIVE) > 0 or len(BENCHMARKS) > 0:
        while len(ACTIVE) < PARALLEL and len(BENCHMARKS) > 0:
            b = BENCHMARKS.pop()
            model = "t_lgl"
            phenotype = "apoptosis"
            bench = f"{model}_{phenotype}"
            # input_file = f"models_phenotype/{bench}"
            output_file = f"{OUT_DIR}/{b}_out.txt"
            command_body = SCRIPT + " " + model + " " + phenotype + " " + "3 " + str(b)
            command = TIMEOUT + " " + CUT_OFF + " time -p " + " " + command_body + " > " + output_file + " 2>&1"
            process = Process(target=SPAWN, args=(command,))
            process.start()
            # pid = pyproc2.find(command_body).pid
            ACTIVE.append((process, bench, output_file))
        time.sleep(1) # Sleep 1s
        STILL_ACTIVE = []
        for (process, bench, output_file) in ACTIVE:
            if process.is_alive():
                STILL_ACTIVE.append((process, bench, output_file))
                # mem_usage = float(subprocess.getoutput(f"ps u -p {pid} | awk '{{sum=sum+$6}}; END {{print sum}}'")) / 1024
                # if mem_usage > MAX_MEM.get(pid, -1):
                #     MAX_MEM[pid] = mem_usage
            else:
                PROCESS_RESULT(process, bench, output_file)
        ACTIVE = STILL_ACTIVE


    # At this point, TIMES should be up-to-date, but we
    # still have to process aggregated data:

    AGGREGATION_LIST = sorted(AGGREGATION_LIST)

    # First entry makes sure the line starts at zero
    AGGREGATION.write(str(AGGREGATION_LIST[0]) + ", 0\n")

    for i in range(len(AGGREGATION_LIST)):
        AGGREGATION.write(str(AGGREGATION_LIST[i]) + ", " + str(i+1) + "\n")

    # Enable this line if you have a 1h timeout.
    #AGGREGATION.write("3600, " + str(len(AGGREGATION_LIST)) + "\n")

    TIMES.close()
    AGGREGATION.close()