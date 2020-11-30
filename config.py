'''
Script to generate Makefile with the correct compiler
configuration for a given machine.
'''

from subprocess import *
import argparse
import sys 

# Manage command-line arguments to load compiler option
parser = argparse.ArgumentParser()
parser.add_argument(
    'config', metavar='CONFIG', type=str,
     help='Name of config file to use (e.g., config/aire_gfortran)')
args = parser.parse_args()


# Determine locations
target_dir  = "./"
config_dir  = "config/"
config_path = args.config 

# Load template Makefile and insert compiler configuration
makefile0    = open(config_dir+"Makefile").read()
compile_info = open(config_path).read()
makefile1 = makefile0.replace("<COMPILER_CONFIGURATION>",compile_info)

# Write the new Makefile to the main directory
open(target_dir+"Makefile","w").write(makefile1)

print( "".join(["\n\nMakefile configuration complete for configuration file: ",config_path,"\n"]) )

instructions = '''==== How to run coordinates ====\n

# Compile a test program
make clean 
make coord-static      # makes the static libray libcoordinates.a

# Run a test program
make test_ccsm3 
./libcoordinates/bin/test_ccsm3.x 

'''

print(instructions)



### Old script commands to automatically determine the config_path from host + compiler:

# Determine the current hostname
#proc = Popen("hostname",shell=True,stdout=PIPE,stderr=PIPE)
#host = proc.communicate()[0].strip().decode("utf-8")

#if host == "manager.cei.ucm.local": host = "eolo"
#if host in ["login01","login02"]:   host = "pik"

#print("{}{}_{}".format(config_dir,host,compiler))

