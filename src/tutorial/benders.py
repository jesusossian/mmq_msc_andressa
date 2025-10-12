#!/usr/bin/python
# ---------------------------------------------------------------------------
# File: benders.py
# Version 22.1.1
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2009, 2022. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ---------------------------------------------------------------------------
"""
Read in a model from a file and solve it using Benders decomposition.

If an annotation file is provided, use that annotation file.
Otherwise, auto-decompose the problem and dump the annotation
to the file 'benders.ann'.

To run this example from the command line, use

    python benders.py filename [annofile]
"""
import sys

import cplex
from cplex.exceptions import CplexError

def create_annotation(cpx):
    """Generates a default annotation.

    Setup and install a default benders partition whereby all
    continuous variables go into a single worker and all other
    variables go into the master partition.
    """
    anno = cpx.long_annotations
    idx = anno.add(name=anno.benders_annotation, defval=anno.benders_mastervalue)
    ctypes = cpx.variables.get_types()
    objtype = anno.object_type.variable
    continuous = cpx.variables.type.continuous
    cpx.long_annotations.set_values(idx, objtype,
                                    [(i, anno.benders_mastervalue + 1)
                                     for i, j
                                     in enumerate(ctypes)
                                     if j == continuous])


def benders(filename, annofile):
    try:
        # Create the Cplex object.
        cpx = cplex.Cplex()

        # Read the problem file.
        cpx.read(filename)

        # If provided, read the annotation file.
        if annofile is not None:
            # Generate default annotations if annofile is "create".
            if annofile == "create":
                create_annotation(cpx)
            else:
                # Otherwise, read the annotation file.
                cpx.read_annotations(annofile)
        else:
            # Set benders strategy to auto-generate a decomposition.
            cpx.parameters.benders.strategy.set(cpx.parameters.benders.strategy.values.full)
            # Write out the auto-generated annotation.
            cpx.write_benders_annotation("benders.ann")
            
        # Add time limit
        cpx.parameters.timelimit.set(360)


        start_time = cpx.get_time()
        # Solve the problem using Benders' decomposition.
        cpx.solve()
        end_time = cpx.get_time()

        # Get the solution status.
        solstatval = cpx.solution.get_status()
        solstatstr = cpx.solution.get_status_string()

        # Get the best bound.
        dualbound = cpx.solution.MIP.get_best_objective()

        # Get the objective function value.
        primalbound = cpx.solution.get_objective_value()
        
        if cpx.solution.get_status() == cpx.solution.status.MIP_OPTIMAL or \
        cpx.solution.get_status() == cpx.solution.status.MIP_FEASIBLE:
            mip_gap_rel = cpx.solution.MIP.get_mip_relative_gap()
            print(f"Relative MIP Gap: {mip_gap_rel}")
        else:
            print("MIP solution not found or not optimal.")
            
        best_integer_objective = c.solution.MIP.get_best_objective()
        best_node_objective = c.solution.MIP.get_best_bound()
        absolute_mip_gap = abs(best_integer_objective - best_node_objective)
        print(f"Absolute MIP Gap: {absolute_mip_gap}")

        # Print the results.
        #print("Solution status: {0}: {1}".format(solstatval, solstatstr))
        #print("Best bound: {0}".format(dualbound))
        #print("Best integer: {0}".format(primalbound))
        
        print(f"Solution status: {solstatval} : {solstatstr}")
        print(f"Best bound: {dualbound}")
        print(f"Best integer: {primalbound}")
        print("Total solve time (sec.):", end_time - start_time)

    except CplexError as exc:
        raise


def usage():
    print("""\
Usage: benders.py filename [annofile]
  filename   Name of a file, with .mps, .lp, or .sav
             extension, and a possible, additional .gz
             extension
  annofile   Optional .ann file with model annotations.
             If "create" is used, the annotation is computed.
Exiting...""")


if __name__ == "__main__":
    # Check the arguments.
    argc = len(sys.argv)
    if argc == 2:
        filename = sys.argv[1]
        annofile = None
    elif argc == 3:
        filename = sys.argv[1]
        annofile = sys.argv[2]
    else:
        usage()
        sys.exit(-1)
        
    # Call the benders function with the appropriate arguments.
    benders(filename, annofile)
