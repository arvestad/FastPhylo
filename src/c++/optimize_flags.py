#!/usr/bin/env python
#

"""
This is a script for searching for the optimal combination of optimization flags for gcc.
The program is executed as:

./optimize_flags -- comandline_for_prog --

The script will then recompile with different combinations of optimization flags and execute
the command line to find the best combination.

The compilation:
1) The flags that the script want to test are defined in the evironment variable CFLAGS
2) make clean is called
3) make is called.

The run:
4) the command line is executed and timed.
"""

import sys
import os
import string
import timeit
import heapq


#
# Each compilation and execution is stored in an instance of the Compilation class
#
class Execution:
    """A class containing the information of the execution of one specific compilation"""
    
    def __init__(self,compilation_flags,execution_time):
        self.compilation_flags = compilation_flags
        self.execution_time = execution_time

    def __cmp__(self, other):
        print type(other)
        if self.execution_time == other.execution_time:
            return 0
        elif self.execution_time < other.execution_time:
            return -1
        else:
            return 1

    def __eq__(self,other):
        return self.compilation_flags == other.compilation_flags
    
    def __hash__(self):
        return hash(self.compilation_flags)
    
    def __str__(self):
        return "execution time: " + str(self.execution_time) + "\nflags: " + self.compilation_flags

#####################################################################
# SEARCH FOR THE BEST COMBINATION OF FLAG_GROUPS
#


#
# A set of instances of Execution
# 
executions = dict()

#
#string with the command line for the program
#
command_line = ""

#
# The groups of optimization flags that we are trying to optimize
#
flag_groups = [] #list of compilation strings for each group

def set_up_groups(num_groups):
    """Makes the global flag_groups list contain num_groups strings,
    each representing a set of compilation flags."""
    global flag_groups
    flag_groups = []
    for i in range(num_groups):
        flag_groups.append("")
    
def add_flag(groupid, flag):
    """Adds a flag to the associated group."""
    global flag_groups
    flag_groups[groupid] += " "+flag


def execute_system_string(s):
    """Executes the input string in the shell."""
    print "Executing: ",s
    sys.stdout.flush()
    os.system(s)
    
#
# Depth first traversal of the combination
#
def depth_first_search(current_group_id=0,current_flags=""):
    """A depth first traversal of all combinations of flags.
    The result is entered into the global dictionary executions."""

    global flag_groups
    global command_line
    #if all groups have been consider compile and execute
    if current_group_id == len(flag_groups):
        print "---------------------------------"
        execution = Execution(current_flags,-1)

        if execution in executions:
            print "compilation allready checked...."
            print executions[execution]
            return
        
        print "setting environment variable CFLAGS = ",current_flags
        os.putenv("CFLAGS",current_flags)
        sys.stdout.flush()
        execute_system_string("rm -f tmp_XXX_dummymake")
        execute_system_string("make clean &>tmp_XXX_dummymake")
        execute_system_string("make &>tmp_XXX_dummymake")
        print "Executing command line: ", command_line
        t = timeit.Timer("os.system(\""+command_line+"\")","import os")       # outside the try/except
        try:        
            sys.stdout.flush()
            execution.execution_time = t.timeit(1)
            executions[execution] = execution
            print "execution time: ", execution.execution_time
            sys.stdout.flush()
        except:
            t.print_exc()
    #try both adding and not adding the flags of the current group
    else:
        depth_first_search(current_group_id+1,current_flags+" "+ flag_groups[current_group_id])
        depth_first_search(current_group_id+1,current_flags)
      


#
# Find the optimal combination of the optimization flags included in "-O1" and "-O2"
#
def first_stage_optimization_flags():
    """Searches for the best combination of flags using the -O1 flag
    and all flags inlcuded in the -O2 option. Each compilation is stored
    in the global variable executions."""
    set_up_groups(5)
    #set up the flag groups    
    add_flag(0, "-O1")
    #02 flags
    add_flag(2,"-fforce-mem")
    add_flag(2,"-foptimize-sibling-calls")
    add_flag(2,"-fstrength-reduce")
    add_flag(3,"-fcse-follow-jumps  -fcse-skip-blocks")
    add_flag(3,"-frerun-cse-after-loop  -frerun-loop-opt")
    add_flag(3,"-fgcse   -fgcse-lm   -fgcse-sm")
    add_flag(2,"-fdelete-null-pointer-checks")
    add_flag(2,"-fexpensive-optimizations")
    add_flag(2,"-fregmove")
    add_flag(4,"-fschedule-insns  -fschedule-insns2")
    add_flag(4,"-fsched-interblock -fsched-spec")
    add_flag(4,"-fcaller-saves")
    add_flag(2,"-fpeephole2")
    add_flag(2,"-freorder-blocks  -freorder-functions")
    add_flag(2,"-fstrict-aliasing")
    add_flag(1,"-falign-functions  -falign-jumps")
    add_flag(1,"-falign-loops  -falign-labels")

    depth_first_search()

#
# Optimization of flags from "-O3"
#
def second_stage_optimization_groups(basic_group):
    """Finds the best combination of flag groups where the supplied
    input string is one of the flag groups. The rest of the flags are
    those included in the -O3 option."""
    set_up_groups(5)
    
    add_flag(0, basic_group)
    
    #O3 flags
    add_flag(3,"-finline-functions")
    add_flag(3,"-frename-registers")
    #
    add_flag(2,"-fno-default-inline")
    add_flag(1,"-fforce-addr")
    add_flag(1,"-fomit-frame-pointer")
    add_flag(1,"-fno-branch-count-reg")
    add_flag(1,"-fstrength-reduce")
    #these flags cause problems for gcc version 3.4.6 
    #add_flag(1,"-fbranch-probabilities")
    #add_flag(1,"-fno-guess-branch-probability")
    add_flag(2,"-ffast-math")
    #
    add_flag(4,"-ftracer")
    add_flag(4,"-funroll-loops")
    add_flag(4,"-fprefetch-loop-arrays")

    depth_first_search()

#
# Not used at the moment. Modify this function to find other combinations and add it to
# the main functio.
#
def third_stage():
    set_up_groups(6)

    add_flag(0,"-fforce-addr")
    add_flag(1,"-fomit-frame-pointer")
    #add_flag(1,"-fno-branch-count-reg")

    #add_flag(1,"-fstrength-reduce")
    #add_flag(2,"-fno-default-inline")
    #2-ffast-math
    add_flag(2,"-finline-functions")
    add_flag(3,"-frename-registers")
    add_flag(4,"-ftracer")
    #add_flag(5,"-funroll-loops")
    add_flag(5,"-fprefetch-loop-arrays")
 
    depth_first_search()
    
def print_usage():
    print "Usage ./optimize_flags -- command_line --"
    sys.exit(1)
    
def main():
    global command_line
    global executions

    executions = dict()
    
    #find the command line
    starti,endi=-1,-1
    for i in range(len(sys.argv)):
        if sys.argv[i] == "--":
            starti=i
            break
    for i in range(starti+1,len(sys.argv)):
        if sys.argv[i] == "--":
            endi=i
            break
    if endi==-1 or starti==-1:
        print_usage()
        
    command_line = string.join(map(str,sys.argv[(starti+1):endi]))
    print "Command line: ", command_line

    #check all combinations in the first stage.
    first_stage_optimization_flags()

    #Use the best combinations from the first stage and combine them with
    #flag groups from the second stage.
    executions_first_stage= [x for x in executions]
    executions_first_stage.sort()
    for i in range(min([5,len(executions_first_stage)])):
        second_stage_optimization_groups(executions_first_stage[i].compilation_flags)

    #Print all executions in order of execution times.
    for x in ["************************************"]*5: print x
    all_executions = [x for x in executions]
    print all_executions
    all_executions.sort()
    tmprange = range(len(all_executions))
    tmprange.reverse()
    for i in tmprange:
        print "----\n",all_executions[i]
    

# If this script is run as a program:    
if __name__=="__main__":
    main()
