import csv
import glob
import os
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
args = vars(ap.parse_args())

rosetta_mini_app = "/work/05388/kiruba/stampede2/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease"
rosetta_db = "/work/05388/kiruba/stampede2/Rosetta/main/database"
rosetta_mini_score_app = "/work/05388/kiruba/stampede2/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease"

# rosetta_mini_app = "/home/kiruba/softwares/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease"
# rosetta_db = "/home/kiruba/softwares/Rosetta/main/database"
# rosetta_mini_score_app = "/home/kiruba/softwares/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease"

# Second minimization of protein ligand complex

def second_minimization(rosetta_mini_app, rosetta_db, protein_ligand_path, lig_params_path):
    dict_prtn_lig_with_params = {}
    for p_l_path in protein_ligand_path:
        params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[4:-1])
        params_file_with_ext = "{0}.params".format(params_file)
        dict_prtn_lig_with_params[p_l_path] = params_file_with_ext

    for key, value in dict_prtn_lig_with_params.items():
        protein_name = os.path.basename(key)
        protein_ligand_without_ext = protein_name.split('.')[0]
        for lig_params in lig_params_path:
            lig_params_basename = os.path.basename(lig_params)
            if lig_params_basename == value:
                print('{0} -database {1} -ignore_unrecognized_res -s {2}/{3} -extra_res_fa {2}/{4} > {2}/mini_protein_ligand_complex_of_all_models/mini_{5}.log'.format(rosetta_mini_app, rosetta_db, os.getcwd(), key, lig_params, protein_ligand_without_ext))
#                os.system('{0} -database {1} -ignore_unrecognized_res -s {2}/{3} -extra_res_fa {2}/{4} > {2}/mini_protein_ligand_complex_of_all_models/mini_{5}.log'.format(rosetta_mini_app, rosetta_db, os.getcwd(), key, lig_params, protein_ligand_without_ext))
        
# Second minimization score_only

def second_minimization_score_only(rosetta_mini_score_app, mini_protein_ligand_path, lig_params_path):
    dict_prtn_lig_with_params = {}
    for p_l_path in mini_protein_ligand_path:
        params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[5:-1])
        params_file_with_ext = "{0}.params".format(params_file)
        dict_prtn_lig_with_params[p_l_path] = params_file_with_ext

    for key, value in dict_prtn_lig_with_params.items():
        protein_name = os.path.basename(key)
        protein_ligand_without_ext = protein_name.split('.')[0]
        for lig_params in lig_params_path:
            lig_params_basename = os.path.basename(lig_params)
            if lig_params_basename == value:
                print('{0} -ignore_unrecognized_res -s {1}/{2} -extra_res_fa {1}/{3} > {1}/mini_protein_ligand_complex_of_all_models/{4}.score'.format(rosetta_mini_score_app, os.getcwd(), key, lig_params, protein_ligand_without_ext))
#                os.system('{0} -ignore_unrecognized_res -s {1}/{2} -extra_res_fa {1}/{3} > {1}/mini_protein_ligand_complex_of_all_models/{4}.score'.format(rosetta_mini_score_app, os.getcwd(), key, lig_params, protein_ligand_without_ext))
        
currentWD = os.getcwd()
os.chdir(args['folder'])

#os.mkdir("mini_protein_ligand_complex_of_all_models")
#second_minimization(rosetta_mini_app, rosetta_db, glob.glob('protein_ligand_complex_top_10_comp_models/*.pdb'), glob.glob('lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.params'))
#os.system("for i in mini_*.pdb; do mv $i; mini_protein_ligand_complex_of_all_models; done")

second_minimization_score_only(rosetta_mini_score_app, glob.glob('mini_protein_ligand_complex_of_all_models/*.pdb'), glob.glob('lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.params'))

os.chdir(currentWD)
