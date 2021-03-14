#!/home/kiruba/softwares/miniconda3/bin/python

import os
import glob
import pandas as pd
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
ap.add_argument("-omega", "--omega_path", required=True,
                help="Enter the OMEGA path")
ap.add_argument("-rocs", "--rocs_path", required=True,
                help="Enter the ROCS path")
ap.add_argument("-temp_lig", "--template_ligand_path", required=True,
                help="Enter the template ligand PDB path")
ap.add_argument("-mol2params", "--mol2params_path", required=True,
                help="Enter the Rosetta mol2params file path")
ap.add_argument("-convert", "--convert_path", required=True,
                help="Convert SDF to mol2 using openeye script")
args = vars(ap.parse_args())


OMEGA = args["omega_path"]
ROCS = args["rocs_path"]
template_lig_library = glob.glob(f"{args['template_ligand_path']}/*.pdb")
mol2params = f"python {args['mol2params_path']}"
convert_py = args["convert_path"]

# Conformer generation using OMEGA openeye (added new flag for resolving the missing MMFF parameters in ligand using -strictatomtyping false)


def conf_gen(smi, maxconfs=2000):

    smi_prefix = os.path.splitext(os.path.basename(smi))[0]
    print('{0} -in {1}/{2} -out {1}/OMEGA/{3}_omega.sdf -prefix {1}/OMEGA/{3}_omega -warts true -maxconfs {4} -strict false'.format(
        OMEGA, os.getcwd(), smi, smi_prefix, maxconfs))
    os.system('{0} -in {1}/{2} -out {1}/OMEGA/{3}_omega.sdf -prefix {1}/OMEGA/{3}_omega -warts true -maxconfs {4} -strict false'.format(
        OMEGA, os.getcwd(), smi, smi_prefix, maxconfs))

# ligand alignment using ROCS openeye


def lig_alignment(conformer, template_database, rocs_maxconfs_output=100):

    sdf_prefix = os.path.basename(os.path.splitext(conformer)[0]).split('_')[0]
    for template in template_database:
        template_id = "_".join(os.path.basename(template).split("_")[0:3])
        # print('{0} -dbase {1}/{2} -query {3} -prefix {4}_{5}_rocs -oformat sdf -maxconfs 30 -outputquery false -qconflabel title -outputdir {1}/ROCS/'.format(ROCS, os.getcwd(),conformer, template, sdf_prefix, template_id))
        os.system('{0} -dbase {1}/{2} -query {3} -prefix {4}_{5}_rocs -oformat sdf -maxconfs 30 -outputquery false -qconflabel title -outputdir {1}/ROCS/'.format(
            ROCS, os.getcwd(), conformer, template, sdf_prefix, template_id))

# Combine each report file into single file


def combine_report_files(report_file):

    dataframeList = []
    for report in report_file:
        target_template_file_name = os.path.basename(
            report).replace("_1.rpt", "")
        read_rpt = pd.read_csv(report, sep='\t', dtype=str)
        read_rpt = read_rpt.loc[:, ~read_rpt.columns.str.match('Unnamed')]
        shapequery_changed_table = read_rpt.replace(
            'untitled-query-1', target_template_file_name)
        shapequery_changed_table.to_csv(report, sep='\t', index=None)
        dataframeList.append(shapequery_changed_table)
    single_rpt_csv_file = pd.concat(dataframeList, axis=0)
    tanimoto_sorted_table = single_rpt_csv_file.sort_values(
        by=['TanimotoCombo'], ascending=False)
    tanimoto_sorted_table.to_csv(
        'ROCS/single_report_file_sorted.csv', index=False)
    top_100_hit_table = tanimoto_sorted_table.iloc[:100, :]
    top_100_conf_temp_names = top_100_hit_table['Name'].astype(
        str) + '_' + top_100_hit_table['ShapeQuery'].astype(str)
    top_100_conf_temp_names.to_csv('ROCS/top_100.txt', sep='\t', index=False)

# Seperate the top 100 conformer hits from ROCS alignment sdf files


def sep_hits_from_rocs_sdf_file(top_100_hits_txt_path):

    with open(top_100_hits_txt_path) as open_file:
        dict_sdf_file_name = {}
        for line in open_file:
            if len(line.split('_')) == 8:
                conf_name = "_".join(line.strip().split('_')[0:3])
                temp_name = "_".join(line.strip().split('_')[3:])
                if conf_name not in dict_sdf_file_name:
                    dict_sdf_file_name[conf_name] = [temp_name]
                else:
                    dict_sdf_file_name[conf_name].append(temp_name)
            elif len(line.split('_')) == 7:
                conf_name = "_".join(line.strip().split('_')[0:2])
                temp_name = "_".join(line.strip().split('_')[2:])
                if conf_name not in dict_sdf_file_name:
                    dict_sdf_file_name[conf_name] = [temp_name]
                else:
                    dict_sdf_file_name[conf_name].append(temp_name)
        for key, list_value in dict_sdf_file_name.items():
            for single_value in list_value:
                print(
                    'sed -n /{0}$/,/\$\$\$\$/p ROCS/{1}_hits_1.sdf > top_100_conf/{0}_{1}_hits.sdf'.format(key, single_value))
                os.system(
                    'sed -n /{0}$/,/\$\$\$\$/p ROCS/{1}_hits_1.sdf > top_100_conf/{0}_{1}_hits.sdf'.format(key, single_value))

# Convert SDF file to PDB/PARAMS for Rosetta input


def sdftomol2(mol2params, top_hits_sdf_path, convert_py):

    for file in top_hits_sdf_path:
        out_put_file_name = os.path.splitext(os.path.basename(file))[0]
        print('python {0} {1} {2}/top_100_conf/{3}.mol2'.format(convert_py, file, os.getcwd(), out_put_file_name))
        os.system('python {0} {1} {2}/top_100_conf/{3}.mol2'.format(convert_py, file, os.getcwd(), out_put_file_name))

        print('{0} -s {1}.mol2 --prefix=mol2params/{2}'.format(mol2params,
                                                               file.split(".")[0], out_put_file_name))
        os.system('{0} -s {1}.mol2 --prefix=mol2params/{2}'.format(mol2params,
                                                                   file.split(".")[0], out_put_file_name))


currentWD = os.getcwd()
os.chdir(args['folder'])

os.mkdir("OMEGA")
conf_gen(glob.glob("*.smi")[0], maxconfs=1000)

os.mkdir("ROCS")
lig_alignment(glob.glob("OMEGA/*.sdf")
              [0], template_lig_library, rocs_maxconfs_output=30)

combine_report_files(glob.glob("ROCS/*.rpt"))

os.mkdir("top_100_conf")
sep_hits_from_rocs_sdf_file("ROCS/top_100.txt")

os.mkdir("mol2params")
sdftomol2(mol2params, glob.glob("top_100_conf/*.sdf"), convert_py)
if not os.listdir("mol2params"):
    os.system("mv *hits_0001.pdb mol2params")
    os.system("mv *hits.params mol2params")
os.chdir(currentWD)
