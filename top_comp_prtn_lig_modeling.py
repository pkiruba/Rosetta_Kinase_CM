import os
import glob
import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
args = vars(ap.parse_args())

# Coping the top 10 ligand pdb and params file

def copy_lig_from_sdf2params(top_10_comp_models, sdf_to_params_path, top_10_mini_models_sdf_to_params):

    for path in top_10_comp_models:
        lig_pdb = '_'.join(os.path.basename(path).split('_')[5:])
        lig_params = '{0}.params'.format('_'.join(os.path.basename(path).replace('_0001.pdb', '').split('_')[5:]))
        print('cp {0}/{1} {2}'.format(sdf_to_params_path, lig_pdb, top_10_mini_models_sdf_to_params))
        os.system('cp {0}/{1} {2}'.format(sdf_to_params_path, lig_pdb, top_10_mini_models_sdf_to_params))
        print('cp {0}/{1} {2}'.format(sdf_to_params_path, lig_params, top_10_mini_models_sdf_to_params))
        os.system('cp {0}/{1} {2}'.format(sdf_to_params_path, lig_params, top_10_mini_models_sdf_to_params))

# Protein Ligand concatenation of all the models

def protein_ligand_concatenation_top_10_comp_models(template_hits, target_fasta_file, ligands):
    
    top_10_comp_models = []
    with open(template_hits) as read_CSV:
        reader = csv.DictReader(read_CSV)
        for row in reader:
            target_pdb_path = os.path.dirname(template_hits)
            target_pdb_name = os.path.basename(target_fasta_file).split('.')[0]
            model = '{0}_{1}.pdb'.format(target_pdb_name, row['template'])
            top_10_comp_models.append(model)

    for model in top_10_comp_models:
        for lig in ligands:
            protein_model = os.path.join(target_pdb_path, model)
            basename_ligand = os.path.splitext(os.path.basename(lig))[0]
            complex_protein_ligand = '{0}_{1}.pdb'.format(model.split('.')[0], basename_ligand)
            print('cat {0} {1} > protein_ligand_complex_top_10_comp_models/{2}'.format(protein_model, lig, complex_protein_ligand))
            os.system('cat {0} {1} > protein_ligand_complex_top_10_comp_models/{2}'.format(protein_model, lig, complex_protein_ligand))

currentWD = os.getcwd()
os.chdir(args['folder'])

os.mkdir("lig_from_top_10_min_models_from_first_comp_model_sdf_to_params")
copy_lig_from_sdf2params(glob.glob('top_10_mini_models_from_first_comp_model/*.pdb'), 'sdf2params', 'lig_from_top_10_min_models_from_first_comp_model_sdf_to_params')

os.mkdir("protein_ligand_complex_top_10_comp_models")            
protein_ligand_concatenation_top_10_comp_models("protein_comp_modeling/template_hits.csv", glob.glob("*.fasta")[0], glob.glob("lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.pdb"))

os.chdir(currentWD)
