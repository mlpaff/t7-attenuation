# Run Kallisto on multiple files within a directory

import argparse
import subprocess
import shlex
import os

def args():
  # define command line options
  parser = argparse.ArgumentParser(description = "Kallisto parameters")
  parser.add_argument('-i', '--index', 
    help="Path to kallisto index", 
    required=True)
  parser.add_argument('-d', '--directory', 
    help='Input directory', 
    required=True)
  parser.add_argument('-l', '--length',
    help='Estimated average fragment length of reads',
    required=True)
  parser.add_argument('-s', '--std_dev',
    help="Estimated standard devition of fragment length",
    required=True)

  args = parser.parse_args()

  return(args)


def main():
  arg = args()
  idx = arg.index
  os.chdir(arg.directory)
  files = [x for x in os.listdir() if x.endswith('.gz')]

  call = "kallisto quant -i "

  for l in range(len(files)):
    cmd = call + idx + ' -o ' + files[l][:-9] + ' --single -l ' + arg.length + ' -s ' + arg.std_dev + ' ' + files[l]
    subprocess.call(cmd, shell=True)
  
main()