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
from shlex import quote

# Spawn a new external process.
def SPAWN(command):
    print("Spawn:", command)
    result = os.system(command)
    sys.exit(result)

# Regex used for parsing the standardised output of
# UNIX time utility (when invoked with -p argument)
RE_TIME = re.compile("\\s*real\\s*(\\d+\\.?\\d*)\\s*")

if __name__ == "__main__":
    print(">>>>>>>>>> START BENCHMARK RUN")

    CUT_OFF = sys.argv[1]
    print("Timeout:", CUT_OFF)
    BENCH_DIR = sys.argv[2]
    print("Benchmark group directory:", BENCH_DIR)
    SCRIPT = sys.argv[3]
    print("Script:", SCRIPT)
    INTERACTIVE = False
    PARALLEL = 0
    if len(sys.argv) > 4:
        INTERACTIVE = sys.argv[4] == '-i'
        if sys.argv[4] == '-p':
            PARALLEL = int(sys.argv[5])

    print("Interactive:", INTERACTIVE)

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


    # Utility function to check if a given path is a benchmark model.
    def is_bench(benchmark):
        return benchmark.endswith(".aeon")

    # Create output directory
    OUT_DIR = BENCH_DIR.replace("/", "_") + "_" + os.path.basename(SCRIPT)
    if PARALLEL > 0:
        OUT_DIR = OUT_DIR + "_parallel"
    OUT_DIR = "_run_" + OUT_DIR + "_" + str(int(time.time()))
    os.mkdir(OUT_DIR)

    # Create output stats file
    TIMES = open(OUT_DIR + "/" + BENCH_DIR.replace("/", "_") + "_" + os.path.basename(SCRIPT) + "_times.csv", "w")
    TIMES.write("Benchmark, Time[s]\n")

    # Create an aggregated stats file
    AGGREGATION = open(OUT_DIR + "/" + BENCH_DIR.replace("/", "_") + "_" + os.path.basename(SCRIPT) + "_aggregated.csv", "w")
    AGGREGATION.write("Time[s], No. Completed\n")

    # Here, save all runtimes.
    AGGREGATION_LIST = []

    BENCHMARKS = []
    for group in os.listdir(BENCH_DIR):
        if not group.startswith("group"):
            continue
        for file in os.listdir(BENCH_DIR+"/"+group):
            if not file.endswith(".aeon"):
                continue
            BENCHMARKS.append(group+"/"+file)
        
    BENCHMARKS = sorted(BENCHMARKS)

    if SCRIPT.endswith(".py"):
        SCRIPT = "python3 " + SCRIPT

    # Handle data from a finished process. In particular,
    # update AGGREGATION_LIST and TIMES file.
    def PROCESS_RESULT(process, name, output_file):
        print("Finished:", output_file)
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

    if PARALLEL > 0:
        ACTIVE = []
        while len(ACTIVE) > 0 or len(BENCHMARKS) > 0:
            while len(ACTIVE) < PARALLEL and len(BENCHMARKS) > 0:
                bench = BENCHMARKS.pop(0)
                name = os.path.splitext(bench)[0].replace("/", "_")
                input_file = BENCH_DIR + "/" + bench
                output_file = OUT_DIR + "/" + name + "_out.txt"
                command = TIMEOUT + " " + CUT_OFF + " " + SCRIPT + " " + input_file + " > " + output_file + " 2>&1"
                process = Process(target=SPAWN, args=(command,))
                process.start()
                ACTIVE.append((process, name, output_file))
            time.sleep(1) # Sleep 1s
            STILL_ACTIVE = []
            for (process, name, output_file) in ACTIVE:
                if process.is_alive():
                    STILL_ACTIVE.append((process, name, output_file))
                else:
                    PROCESS_RESULT(process, name, output_file)
            ACTIVE = STILL_ACTIVE
    else:
        for bench in BENCHMARKS:
            print(">>>>>>>>>> START BENCHMARK", bench, "<<<<<<<<<<")
            if INTERACTIVE:
                print("Write 'skip' to go to next benchmark, 'abort' to end the run, or press enter key to continue...")
                skip = sys.stdin.readline()
                if skip.startswith("skip"):
                    print("Skipped!")
                    continue
                if skip.startswith("abort"):
                    print("Aborted!")
                    break
            # Benchmark filename without extension
            name = os.path.splitext(bench)[0].replace("/", "_")
            print(name)
            input_file = BENCH_DIR + "/" + bench
            output_file = OUT_DIR + "/" + name + "_out.txt"
            command = TIMEOUT + " " + CUT_OFF + " " + SCRIPT + " " + input_file + " > " + output_file + " 2>&1"
            process = Process(target=SPAWN, args=(command,))
            process.start()
            process.join()
            PROCESS_RESULT(process, name, output_file)


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