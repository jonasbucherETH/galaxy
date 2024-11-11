#!/usr/bin/env python
 
import argparse
import subprocess
import os
import shutil
import zipfile
import re
import tempfile
import sys
from glob import glob

def run_bismark2bedGraph(input_files, output, output_bedgraph, output_coverage, cytosine_context):
    # Build the bismark2bedGraph command

    # Get the directory where the input file is located
    input_file_dir = os.path.dirname(os.path.abspath(input_files))

    prev_dir = os.getcwd()
    output_dir = tempfile.mkdtemp()

    #command = ["bismark2bedGraph", "--dir", output_dir, "--output", output]
    command = ["bismark2bedGraph", "--dir", output_dir, "--output", output]
    #command = ["bismark2bedGraph", "--output", output]

    # Add CX_context flag if enabled
    if cytosine_context=="CpG":
        i = 0
        for file in input_files.split(','):
            new_basename = "CpG" + str(i)
            i+=1
            new_file = os.path.join(output_dir, new_basename)
            shutil.move(file, new_file)
            command.append(new_file)
    else:
        command.append("--CX")
        for file in input_files.split(','):
            command.append(file)
    

    # Output bedGraph file
    #command.extend(["--scaffolds"])

    #shutil.move(glob(os.path.join(output_dir, bed_out_name))[0], output_bedgraph)
    #shutil.move(glob(os.path.join(output_dir, cov_out_name))[0], output_coverage)

    os.chdir(
        output_dir
    )  # needed due to a bug in bismark where the coverage file cannot be found

    #abs_input_files = os.path.abspath(input_files)
    #command.append(os.path.basename(input_file))
    #print("input_files:", input_files)


    # Run the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        raise

#    files_output_dir = os.listdir(output_dir)
#    print("files_output_dir", files_output_dir)
    os.chdir(prev_dir)

    bed_out_name = output + ".gz"
    cov_out_name = output + ".gz.bismark.cov.gz"

#    print("bed_out_name", bed_out_name)
#    print("cov_out_name", cov_out_name)
#    print("output_dir", output_dir)
#    print("prev_dir", prev_dir)
#    print("output_bedgraph", output_bedgraph)
#    print("output_coverage", output_coverage)
#    print("glob(os.path.join(output_dir, bed_out_name)):", glob(os.path.join(output_dir, bed_out_name)))

    #shutil.move(os.path.join(output_dir, bed_out_name), output_bedgraph)
    #shutil.move(os.path.join(output_dir, cov_out_name), output_coverage)
    shutil.copy(glob(os.path.join(output_dir, bed_out_name))[0], output_bedgraph)
    shutil.copy(glob(os.path.join(output_dir, cov_out_name))[0], output_coverage)

#    files_prev_dir = os.listdir(prev_dir)
#    print("files_prev_dir", files_prev_dir)

    #zipper(output_dir, compress)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Galaxy wrapper for bismark2bedGraph')
    parser.add_argument('--input_files', required=True, help='Input file from Bismark methylation extractor')
    parser.add_argument('--output', required=True, help='Output file path name')
    parser.add_argument('--output_bedgraph', required=True, help='Output file bedgraph')
    parser.add_argument('--output_coverage', required=True, help='Output file coverage')
    parser.add_argument('--cytosine_context', required=True, help='Cytosine context')
    
    args = parser.parse_args()

    # Run bismark2bedGraph with the provided arguments
    run_bismark2bedGraph(
        input_files=args.input_files,
        output=args.output,
        output_bedgraph=args.output_bedgraph,
        output_coverage=args.output_coverage,
        cytosine_context=args.cytosine_context
    )
