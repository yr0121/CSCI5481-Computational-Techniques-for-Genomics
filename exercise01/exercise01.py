# Template file for CSCI 5481 Fall 2015 Exercise 1
# Usage:
# template.py -h
# discussed with Huilin Zhang
import sys, os
import argparse
from subprocess import Popen, PIPE

def make_arg_parser():
    parser = argparse.ArgumentParser(prog='template.py',
                          version="%prog 1.0",
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]") 
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-t","--taxonomy",
                      default=None,
                      required=True,
                      help="Path to taxonomy file [required]")
    parser.add_argument("-o","--output",
                      default=None,
                      required=True,
                      help="Path to output file [required]")
    parser.add_argument("-c","--command",
                      default='./burst',
                      help="Path to BURST command") 
    parser.add_argument("-V","--verbose",
                      action="store_true",
                      help="Verbose output")
    return parser

# Runs BURST to search query sequences against reference sequences
def run_burst(query, ref, taxonomy, output, burst_cmd='./burst', verbose=False):
    """thread worker function"""
    #read the command as provided in the homework instruction
    cmd="./burst -r ref.fna -q query.fna --taxonomy taxonomy.txt -o output.txt"
    return run_command(cmd, verbose=verbose)

# runs the given command and returns return value and output
def run_command(cmd, verbose=False):
    if verbose:
        print (cmd)
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return return_value, stdout, stderr

    
if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    #obtain the result from calling the run_burst function and run_command function
    return_value, stdout, stderr=run_burst(args.query,args.ref,args.taxonomy,args.output,args.command,args.verbose)
    #print out the return value based on input -V
    if args.verbose:
        print(return_value)
        print(stdout)
        print(stderr)
    ### YOUR CODE HERE ###

