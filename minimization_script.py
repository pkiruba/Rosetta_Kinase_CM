import os
import glob
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
ap.add_argument("-ma", "--mini_app", required=True, help="Minimization app path")
ap.add_argument("-rf", "--res_file", required=True, help="Residue file path")
args = vars(ap.parse_args())

rosetta_mini_app = args["mini_app"]
residue_of_interest = args["res_file"]
rosetta_mini_score_app = "/work/05388/kiruba/stampede2/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease"

# First minimization of protein ligand complex

def first_minimization(rosetta_mini_app, protein_ligand_path, lig_params_path, residue_of_interest):
    
    dict_prtn_lig_with_params = {}
    for p_l_path in protein_ligand_path:
        params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[4:-1])
        params_file_with_ext = "{0}.params".format(params_file)
        dict_prtn_lig_with_params[params_file_with_ext] = p_l_path

    for params in lig_params_path:
        params_basename = os.path.basename(params)
        if params_basename in dict_prtn_lig_with_params:
            protein_ligand = os.path.basename(dict_prtn_lig_with_params[params_basename])
            protein_ligand_without_ext = os.path.splitext(os.path.basename(dict_prtn_lig_with_params[params_basename]))[0]
            template_id = "_".join(protein_ligand_without_ext.split("_")[2:4])
            print('python {0} -s {1}/{2} -extra_res_fa {1}/{3} -residue_of_interest {4}/{5}_res_file.txt > {1}/mini_protein_ligand_complex_top_1_comp_model/mini_{6}.log'.format(rosetta_mini_app, os.getcwd(), dict_prtn_lig_with_params[params_basename], params, residue_of_interest, template_id, protein_ligand_without_ext))
            

# First minimization score_only

def first_minimization_score_only(rosetta_mini_score_app, mini_protein_ligand_path, lig_params_path):
    
    dict_prtn_lig_with_params = {}
    for mini_p_l_path in mini_protein_ligand_path:
        params_file = "_".join(os.path.splitext(os.path.basename(mini_p_l_path))[0].split("_")[5:-1])
        params_file_with_ext = "{0}.params".format(params_file)
        dict_prtn_lig_with_params[params_file_with_ext] = mini_p_l_path
        
    for params in lig_params_path:
        params_basename = os.path.basename(params)
        if params_basename in dict_prtn_lig_with_params:
            protein_ligand = os.path.basename(dict_prtn_lig_with_params[params_basename])
            protein_ligand_without_ext = os.path.splitext(os.path.basename(dict_prtn_lig_with_params[params_basename]))[0]
            print('{0} -ignore_unrecognized_res -s {1}/{2} -extra_res_fa {1}/{3} > {1}/mini_protein_ligand_complex_top_1_comp_model/{4}.score'.format(rosetta_mini_score_app, os.getcwd(), dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))
            # os.system('{0} -ignore_unrecognized_res -s {1}/{2} -extra_res_fa {1}/{3} > {1}/mini_protein_ligand_complex_top_1_comp_model/{4}.score'.format(rosetta_mini_score_app, os.getcwd(), dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))

currentWD = os.getcwd()
os.chdir(args['folder'])

os.mkdir("mini_protein_ligand_complex_top_1_comp_model")
first_minimization(rosetta_mini_app, glob.glob('protein_ligand_complex_top_1_comp_model/*.pdb'), glob.glob('mol2params/*.params'), residue_of_interest)
#os.system("for i in mini_*.pdb; do mv $i; mini_protein_ligand_complex_top_1_comp_model; done")

# first_minimization_score_only(rosetta_mini_score_app, glob.glob('mini_protein_ligand_complex_top_1_comp_model/*.pdb'), glob.glob('sdf2params/*.params'))

os.chdir(currentWD)
