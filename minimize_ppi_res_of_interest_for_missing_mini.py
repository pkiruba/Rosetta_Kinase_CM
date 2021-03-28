#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license. The Rosetta software is developed by the contributing
# (c) members of the Rosetta Commons. For more information, see
# (c) http://www.rosettacommons.org. Questions about this can be addressed to
# (c) University of Washington CoMotion, email: license@uw.edu.


"""Script of minimization protein-ligand complex. Code created based on the
   minimize_ppi of C++ Rosetta. Original authors: jk, ragul, skelow, bazzoli.
   I added two new arguments - minimization convergence and max numbet of
   iteration. Also I removed some non-important for minimization functions like
   `print unbound`, etc.

   Usage: ./minimize_ppi.py -s <complex> -extra_res_fa <ligand_params_file>

   Author: Grigorii Andrianov
   Small contribution: Kirubakaran Palani (Added the code for minimizing only the Kinase ATP active site residues to reduce the minimization process)
"""

import argparse
import math
import time
import os

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.conformation import Conformation
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.optimization import AtomTreeMinimizer
from pyrosetta.rosetta.core.optimization import MinimizerOptions
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.pose.metrics import CalculatorFactory
from pyrosetta.rosetta.core.pose.metrics.simple_calculators import \
    SasaCalculatorLegacy
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.scoring.constraints import ConstraintSet
from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.protocols.pose_metric_calculators import \
    PackstatCalculator
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import \
    BuriedUnsatisfiedPolarsCalculator as BurUnsatCalc
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import \
    NumberHBondsCalculator


def minimization(bound_pose, scorefxn, convergence, max_iter, pdb_res_num_to_pose_res_num):
    """Minimize input structures using atom tree minimizer

    Args:
        bound_pose (pose): Input pose for minimization
        scorefxn (scorexn): Score function
        convergence (float): Convergence threshold of minimization
        max_iter (int): Max number of minimizations if convergense isn't
                        reachable
        pdb_res_num_to_pose_res_num (list): Selected residues for minimization                

    Returns:
        pose: Return minimized input pose
    """
    minimizer = AtomTreeMinimizer()
    min_options = MinimizerOptions("dfpmin", convergence, True, False)
    min_options.max_iter(max_iter)
    
    if len(pdb_res_num_to_pose_res_num) > 0:
        movemap = MoveMap()
        movemap.set_chi(False)
        movemap.set_bb(False)
        movemap.set_jump(True)

        for res_num in pdb_res_num_to_pose_res_num:
            movemap.set_chi(res_num, True)
            movemap.set_bb(res_num, True)

        minimizer.run(bound_pose, movemap, scorefxn, min_options)
        scorefxn(bound_pose)
        total = bound_pose.energies().total_energies()
        print("Post minimization 1 constrained score: ",
              total[ScoreType.total_score])

        return bound_pose

    else:    
        mm_all = MoveMap()
        mm_all.set_chi(True)
        mm_all.set_bb(True)
        mm_all.set_jump(True)
        minimizer.run(bound_pose, mm_all, scorefxn, min_options)

        

        scorefxn(bound_pose)
        total = bound_pose.energies().total_energies()
        print("Post minimization 1 constrained score: ",
              total[ScoreType.total_score])

        return bound_pose


def calculators(bound_pose, unbound_pose):
    """Printing difference between bound and unbound poses

    Args:
        bound_pose (pose): Bound with ligand pose
        unbound_pose (pose): Unbound with ligand pose (i.e. ligand pushed away
                             from protein)
    """
    sasa_calc_name = "sasa"
    hbond_calc_name = "hbond"
    packstat_calc_name = "packstat"
    burunsat_calc_name = "burunsat"

    calculator_factory = CalculatorFactory.Instance()

    sasa_calculator = SasaCalculatorLegacy()
    calculator_factory.register_calculator(sasa_calc_name, sasa_calculator)

    hb_calc = NumberHBondsCalculator()
    calculator_factory.register_calculator(hbond_calc_name, hb_calc)

    packstat_calc = PackstatCalculator()
    calculator_factory.register_calculator(packstat_calc_name, packstat_calc)

    burunsat_calc = BurUnsatCalc(sasa_calc_name, hbond_calc_name)
    calculator_factory.register_calculator(burunsat_calc_name, burunsat_calc)

    # calculate interface Energy
    bound_energy = bound_pose.energies().total_energies()[
        ScoreType.total_score]
    unbound_energy = unbound_pose.energies().total_energies()[
        ScoreType.total_score]
    interface_energy = bound_energy - unbound_energy

    # delta sasa calculation
    bound_sasa = float(bound_pose.print_metric(
        sasa_calc_name, "total_sasa"))
    unbound_sasa = float(unbound_pose.print_metric(
        sasa_calc_name, "total_sasa"))
    total_bsa = unbound_sasa - bound_sasa

    # interface hb calculation

    bound_hb = float(bound_pose.print_metric(
        hbond_calc_name, "all_Hbonds"))
    unbound_hb = float(unbound_pose.print_metric(
        hbond_calc_name, "all_Hbonds"))
    interface_hb = bound_hb - unbound_hb

    # packstat calculation
    bound_packstat = float(bound_pose.print_metric(
        packstat_calc_name, "total_packstat"))
    unbound_packstat = float(unbound_pose.print_metric(
        packstat_calc_name, "total_packstat"))
    total_packstats = bound_packstat - unbound_packstat

    # unsat polar calculation
    bound_unsat = float(bound_pose.print_metric(
        burunsat_calc_name, "all_bur_unsat_polars"))
    unbound_unsat = float(unbound_pose.print_metric(
        burunsat_calc_name, "all_bur_unsat_polars"))
    interface_unsat = bound_unsat - unbound_unsat

    # Description of output scores
    print("Description of output scores")
    print("Interface_Scores: bound_energy\tinterface_energy\ttotal_bsa\t",
          "interface_hb\ttotal_packstats\tinterface_unsat")

    print("Interface_Scores:", bound_energy, interface_energy, total_bsa,
          interface_hb, total_packstats, interface_unsat, sep="\t")


def ligand_rmsd(pose1, pose2):
    """Calculate RMSD of two ligands

    Args:
        pose1 (pose): Pose 1
        pose2 (pose): Pose 2

    Returns:
        float: RMSD score
    """
    # define ligand residue id of pose 1
    ref_res_num = 0
    for j in range(1, pose1.size() + 1):
        if not pose1.residue(j).is_protein():
            ref_res_num = j
            break

    pose1_rsd = pose1.conformation().residue(ref_res_num)

    inp_res_num = 0

    for j in range(1, pose2.size() + 1):
        if not pose2.residue(j).is_protein():
            inp_res_num = j
            break

    pose2_rsd = pose2.conformation().residue(inp_res_num)

    j = 1
    dist_sum = 0
    for i in range(1, pose1_rsd.nheavyatoms() + 1):
        if j <= pose2_rsd.nheavyatoms():
            x_dist = (pose1_rsd.atom(i).xyz()[
                      0] - pose2_rsd.atom(i).xyz()[0]) ** 2
            y_dist = (pose1_rsd.atom(i).xyz()[
                      1] - pose2_rsd.atom(i).xyz()[1]) ** 2
            z_dist = (pose1_rsd.atom(i).xyz()[
                      2] - pose2_rsd.atom(i).xyz()[2]) ** 2
            dist_sum += x_dist + y_dist + z_dist
        j += 1

    rmsd = math.sqrt(dist_sum / pose1_rsd.nheavyatoms())

    return rmsd


def get_unbound(bound_pose, scorefxn):
    """Push away ligand from protein

    Args:
        bound_pose (pose): Bounded with ligand pose
        scorefxn (scorefxn): Score function

    Returns:
        pose: Unbound pose
    """
    unbound_pose = bound_pose.clone()
    unbound_dist = 80.0

    # use the LAST jump as the one between partners
    rb_jump = bound_pose.num_jump()
    trans_mover = RigidBodyTransMover(unbound_pose, rb_jump)
    trans_mover.trans_axis(trans_mover.trans_axis())
    trans_mover.step_size(unbound_dist)
    trans_mover.apply(unbound_pose)

    scorefxn(unbound_pose)

    return unbound_pose


def main():
    """Main function of minimization
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", type=str, required=True)
    parser.add_argument("-extra_res_fa", type=str, required=True)
    parser.add_argument("-residue_of_interest")
    parser.add_argument("-mini_output", type=str, default=None)
    parser.add_argument("-convergence", type=float, default=0.00001)
    parser.add_argument("-max_iter", type=int, default=2000)
    parser.add_argument("-cst_force_constant", type=float, default=0.0000001)
    parser.add_argument("-debug", type=int, default=None)

    args = parser.parse_args()

    s = args.s
    extra_res_fa = args.extra_res_fa
    residue_of_interest = args.residue_of_interest
    mini_output = args.mini_output
    convergence = args.convergence
    max_iter = args.max_iter
    cst_force_constant = args.cst_force_constant
    debug = args.debug


    # pyrosetta init
    init = " -gen_potential -extra_res_fa {0}".format(extra_res_fa)
    if debug is not None:
        init += " " + "-no_optH true -flip_HNQ false -constant_seed -run:jran \
        {0}".format(debug)

    print(init)
    pyrosetta.init(init)

    # setup scorefxn
    scorefxn = get_score_function()

    # setup the bound pose
    bound_pose = pose_from_pdb(s)
    nres = bound_pose.size()

    pre_min_darc_pose = Pose().assign(bound_pose)

    # initial score
    scorefxn(bound_pose)
    total = bound_pose.energies().total_energies()
    print("Initial score: ", total[ScoreType.total_score])

    # constraints for minimization
    cst_weight = 1
    # take reciprocal and sqrt to pass as force constant
    coord_sdev = math.sqrt(1 / cst_force_constant)

    cst_set = ConstraintSet()
    spring = HarmonicFunc(0.0, coord_sdev)
    conformation = Conformation(bound_pose.conformation())

    my_anchor = 1
    fold_tree = bound_pose.fold_tree()
    rerooted_fold_tree = fold_tree
    rerooted_fold_tree.reorder(my_anchor)
    bound_pose.fold_tree(rerooted_fold_tree)

    for i in range(1, nres + 1):
        rsd = bound_pose.residue(i)
        for j in range(1, rsd.nheavyatoms() + 1):
            atom = AtomID(j, i)
            coord_cst = CoordinateConstraint(atom, AtomID(1, my_anchor),
                                             conformation.xyz(atom), spring)
            cst_set.add_constraint(coord_cst)

    bound_pose.constraint_set(cst_set)
    scorefxn.set_weight(ScoreType.coordinate_constraint, cst_weight)

    start_time = time.time()
    print("Starting minimization....")

    if mini_output is None:
        mini_output = "mini_" + s.split("/")[-1]
    
    pdb_res_num_to_pose_res_num = []
    if residue_of_interest is not None:
        chain_id = os.path.basename(residue_of_interest).split("_")[1]
        with open(residue_of_interest) as f_obj:
            for pdb_res_num in f_obj:
                pdb_res_num = int(pdb_res_num.strip())
                # Here the reason for making this condition is there are some models which has HETATM in template and when 
                # protein modeling part pyrosetta did not considered those residues and in my residue number list I some how keep those
                # residue numbers which results in to 0 when pdb2pose and stops minimization. Now, I look each residue number and 
                # if the residue number is 0 then I am not including that in my minimization part.
                if bound_pose.pdb_info().pdb2pose(chain_id, pdb_res_num) == 0:
                    continue
                else:
                    pdb_res_num_to_pose_res_num.append(bound_pose.pdb_info().pdb2pose(chain_id, pdb_res_num))  

    bound_pose = minimization(bound_pose, scorefxn, convergence, max_iter, pdb_res_num_to_pose_res_num)

    # final minimization
    bound_pose.remove_constraints()

    scorefxn(bound_pose)

    # align minimized pose to the original docked pose and dump pdb complex
    # and ligand
    sp_mover = SuperimposeMover()
    sp_mover.set_reference_pose(pre_min_darc_pose, 1,
                                pre_min_darc_pose.size() - 1)
    sp_mover.set_target_range(1, (bound_pose.size() - 1))
    sp_mover.apply(bound_pose)
    scorefxn(bound_pose)
    bound_pose.dump_scored_pdb(mini_output, scorefxn)

    print("Successfully finished minimizing input.")
    
    scorefxn(bound_pose)
    total = bound_pose.energies().total_energies()
    print("Final score: ", total[ScoreType.total_score])

    # setup the unbound pose
    unbound_pose = get_unbound(bound_pose, scorefxn)

    # print RMSD between input and minimized ligands (only after
    # superposition)
    rmsd = ligand_rmsd(pre_min_darc_pose, bound_pose)
    print("computing post minization ligand RMSD: ", rmsd)

    # define containers for metrics for total complex
    calculators(bound_pose, unbound_pose)

    end_time = time.time()
    print("Total time taken (in seconds): ", round(end_time-start_time, 2))

if __name__ == "__main__":
    main()
