import os
import glob
import pandas as pd
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
args = vars(ap.parse_args())

# First minimization analysis

def analysis_1(score_file_paths):
    
    table_df = pd.DataFrame(columns = ('MODEL', 'BOUND', 'UNBOUND', 'DIFFERENCE'))
    ind=1
    for path in score_file_paths:
        model = os.path.basename(path.split('.')[0])
        with open(path, 'r') as rf:
            for line in rf:
                if 'apps.pilot.bou-min-ubo-nrg-jump: bound pose\'s energy:' in line:
                    bound_pose_energy = line.split()[-1]
                elif 'apps.pilot.bou-min-ubo-nrg-jump: unbound pose\'s energy:' in line:
                    unbound_pose_energy = line.split()[-1]
                elif 'apps.pilot.bou-min-ubo-nrg-jump: energy difference:' in line:
                    energy_difference = line.split()[-1]
                    table_df.loc[ind] = model, float(bound_pose_energy), float(unbound_pose_energy), float(energy_difference)
        ind+=1

    table_df = table_df.sort_values(by =['DIFFERENCE'], ascending=True)
    table_df.to_csv('mini_protein_ligand_complex_top_1_comp_model/rosetta_energy_table.csv', sep=',', index=None)
    table_df['MODEL'][0:10].to_csv('mini_protein_ligand_complex_top_1_comp_model/top_10_models_based_on_rosetta_energy.txt', index=None)
    os.mkdir('top_10_mini_models_from_first_comp_model')
    table_df[0:10].to_csv('top_10_mini_models_from_first_comp_model/top_10_models_rosetta_energy.csv', sep=',', index=None)
    for model in table_df['MODEL'][0:10]:
        pre_mini_model = '{}.pdb'.format('_'.join(model.split('_')[1:]))
        mini_model = '{}.pdb'.format(model)
        print('cp mini_protein_ligand_complex_top_1_comp_model/{} top_10_mini_models_from_first_comp_model'.format(mini_model))
        os.system('cp mini_protein_ligand_complex_top_1_comp_model/{} top_10_mini_models_from_first_comp_model'.format(mini_model))

currentWD = os.getcwd()
os.chdir(args['folder'])

score_file = glob.glob('mini_protein_ligand_complex_top_1_comp_model/*.score')
analysis_1(score_file)

os.chdir(currentWD)
