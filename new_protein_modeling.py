import argparse
import os
import glob
import pandas as pd
from collections import defaultdict
import pyrosetta
pyrosetta.init()

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--folder", required=True, help="Enter into the folder")
ap.add_argument("-emboss", "--embosss_path", required=True,
                help="Enter the EMBOSS path")
ap.add_argument("-temp_seq", "--template_seq_path", required=True,
                help="Enter the template kinase sequence path")
ap.add_argument("-apo_pdb", "--apo_pdb_path", required=True,
                help="Enter the kinase APO PDB path")
args = vars(ap.parse_args())

emboss_needle = args["embosss_path"]
tempate_seqs = glob.glob(f"{args['template_seq_path']}/*.fasta")
apo_pdbs = args["apo_pdb_path"]

# Sequence alignment using EMBOSS


def emboss_needle_search(emboss_needle, target_seq_path, template_seq_path):

    for template_seq in template_seq_path:
        target_seq_id = os.path.basename(target_seq_path).split('.')[0]
        template_seq_id = os.path.basename(template_seq).split('.')[0]
        print('{0} -sid1 {1} -asequence {2}/{3} -sid2 {4} -bsequence {5} -gapopen 10.0 -gapextend 0.5 -aformat3 markx3 -outfile {2}/protein_comp_modeling/protein_seq_alignment_files/{1}_{4}.needle'.format(
            emboss_needle, target_seq_id, os.getcwd(), target_seq_path, template_seq_id, template_seq))
        os.system('{0} -sid1 {1} -asequence {2}/{3} -sid2 {4} -bsequence {5} -gapopen 10.0 -gapextend 0.5 -aformat3 markx3 -outfile {2}/protein_comp_modeling/protein_seq_alignment_files/{1}_{4}.needle'.format(
            emboss_needle, target_seq_id, os.getcwd(), target_seq_path, template_seq_id, template_seq))

# Select top 10 hit templates based on Template Score method


def select_top_hits_from_emboss_and_rocs_pdb(emboss_align_file_path, rocs_align_file_path, target_seq_path):

    emboss_result = pd.DataFrame(columns=(
        'query', 'template', 'length', 'identity', 'similarity', 'gaps', 'score'))
    emboss_ind = 1
    for emboss_alignment_file in emboss_align_file_path:
        with open(emboss_alignment_file, 'r') as rf:
            for line in rf:
                line = line.strip()
                if '# 1:' in line:
                    query = line.split()[-1]
                elif '# 2:' in line:
                    template = line.split()[-1]
                elif '# Length:' in line:
                    length = line.split()[-1]
                elif '# Identity:' in line:
                    identity = line.split(
                    )[-1].replace('(', '').replace(')', '').replace('%', '')
                elif '# Similarity:' in line:
                    similarity = line.split(
                    )[-1].replace('(', '').replace(')', '').replace('%', '')
                elif '# Gaps:' in line:
                    gaps = line[24:].replace(
                        '(', '').replace(')', '').replace('%', '')
                elif '# Score:' in line:
                    score = line.split(
                    )[-1].replace('(', '').replace(')', '').replace('%', '')
#             print (query, template, length, identity, similarity, gaps, score)
            emboss_result.loc[emboss_ind] = query, template, float(length), float(
                identity), float(similarity), float(gaps), float(score)
            emboss_ind += 1
    # Here ranking the df in descending (ascending=False) order and if found two scores are same
    # then assign the same rank. Dense rank does not skip any rank (in min and max ranks are skipped))
    emboss_result['score_rank'] = emboss_result['score'].rank(
        method='dense', ascending=False)
    rank_emboss_hits_dict = emboss_result.set_index(
        "template")["score_rank"].to_dict()

    rocs_single_rpt_df = pd.read_csv(rocs_align_file_path, sep=",")
    rocs_single_rpt_df["row_name_and_shapequery_conc"] = rocs_single_rpt_df["Name"].str.cat(
        rocs_single_rpt_df["ShapeQuery"], sep="_")
    rocs_result = rocs_single_rpt_df.loc[:, ["row_name_and_shapequery_conc", "TanimotoCombo"]].rename(
        columns={"row_name_and_shapequery_conc": "rocs_pdb"})
    rocs_result['TanimotoCombo_rank'] = rocs_result['TanimotoCombo'].rank(
        method='dense', ascending=False)
    rank_rocs_hits_dict = rocs_result.set_index(
        "rocs_pdb")["TanimotoCombo_rank"].to_dict()

    final_hit_list_df = pd.DataFrame(columns=(
        'template', 'score_rank', 'rocs_pdb', 'TanimotoCombo_rank', 'final_rank'))
    ind_title = 1
    for template, score_rank in rank_emboss_hits_dict.items():
        for rocs_pdb, TanimotoCombo_rank in rank_rocs_hits_dict.items():
            if len(rocs_pdb.split('_')) == 8:
                pdb_with_chain = '_'.join(rocs_pdb.split('_')[4:6])
                target_seq_pdb = os.path.basename(
                    target_seq_path).split('.')[0]
                if pdb_with_chain == template and target_seq_pdb != template:
                    final_hit_list_df.loc[ind_title] = template, score_rank, pdb_with_chain, TanimotoCombo_rank, score_rank*TanimotoCombo_rank
                    ind_title += 1
            elif len(rocs_pdb.split('_')) == 7:
                pdb_with_chain = '_'.join(rocs_pdb.split('_')[3:5])
                target_seq_pdb = os.path.basename(
                    target_seq_path).split('.')[0]
                if pdb_with_chain == template and target_seq_pdb != template:
                    final_hit_list_df.loc[ind_title] = template, score_rank, pdb_with_chain, TanimotoCombo_rank, score_rank*TanimotoCombo_rank
                    ind_title += 1
    final_hit_list_df = final_hit_list_df.drop_duplicates(
        subset='rocs_pdb', keep="first")
    final_hit_list_df = final_hit_list_df.sort_values(
        'final_rank', ascending=True)
    final_hit_list_df.iloc[:10].to_csv(
        'protein_comp_modeling/template_hits.csv', sep=',', index=False)

# Top 10 Modeling of target protein using 10 templates


def modeling(template_hits, template_pdb_path, alignment_file_path, target_seq_path):

    templates = pd.read_csv(template_hits, sep=",")
    target_seq = os.path.basename(target_seq_path).split('.')[0]
    templates['tar_tem_seq_alin'] = templates['template'].apply(
        lambda x: "{}_{}.needle".format(target_seq, x))
    templates['tar_tem_seq_alin'] = alignment_file_path + \
        templates['tar_tem_seq_alin']
    top_hit_template_file_path = templates['tar_tem_seq_alin'].tolist()

    aligned_seq = defaultdict(list)
    for path in top_hit_template_file_path:
        target_template_file_name = os.path.splitext(os.path.basename(path))[0]
        target_name_fasta_format = '>{} ..'.format(
            '_'.join(target_template_file_name.split('_')[0:2]))
        template_name_fasta_format = '>{} ..'.format(
            '_'.join(target_template_file_name.split('_')[2:]))
        target_aligned_seq = ''
        template_aligned_seq = ''
        with open(path, 'r') as readFile:
            parse = False
            parse2 = False
            for line in readFile:
                line = line.strip()
                if not parse:
                    if line.startswith(target_name_fasta_format):
                        parse = True
                elif line.startswith(template_name_fasta_format):
                    parse = False
                else:
                    target_aligned_seq += line

                if not parse2:
                    if line.startswith(template_name_fasta_format):
                        parse2 = True
                elif line.startswith('#'):
                    parse2 = False
                else:
                    template_aligned_seq += line
        aligned_seq[target_template_file_name].append(target_aligned_seq)
        aligned_seq[target_template_file_name].append(template_aligned_seq)

    target_seq_for_modeling = {}
    for name, alignment_file in aligned_seq.items():
        top_hits_alignment = '{}\n{}\n{}\n\n'.format(
            name, alignment_file[0], alignment_file[1])
        with open('protein_comp_modeling/top_hits_alignment.txt', 'a') as writeFile:
            writeFile.write(top_hits_alignment)
        target_seq_based_on_temp_pdb = ''
        for i in range(len(alignment_file[0])):
            if not alignment_file[1][i] == '-':
                target_seq_based_on_temp_pdb += alignment_file[0][i]
        target_seq_for_modeling[name] = target_seq_based_on_temp_pdb

    final_target_template_for_modeling = {}
    for target_template, target_final_seq in target_seq_for_modeling.items():
        template_name = '_'.join(target_template.split('_')[2:])
        temp_list_dir = os.listdir(template_pdb_path)
        for template_hit in temp_list_dir:
            if template_name in template_hit:
                final_target_template_for_modeling[template_hit] = target_final_seq

    for template_pdb, target_seq in final_target_template_for_modeling.items():
        output_model_name = 'protein_comp_modeling/{}_{}.pdb'.format(
            target_seq_path.split('.')[0], '_'.join(template_pdb.split('_')[0:2]))
        join_apo_dir_path = os.path.join(template_pdb_path, template_pdb)
        pose = pyrosetta.pose_from_file(join_apo_dir_path)
        assert(pose.size() == len(target_seq))
        scorefxn = pyrosetta.get_fa_scorefxn()
        for i in range(len(target_seq)):
            seqpos = i + 1
            name1 = target_seq[i]
            if (name1 == "-"):
                continue
            pyrosetta.rosetta.protocols.toolbox.pose_manipulation.repack_this_residue(
                seqpos, pose, scorefxn, True, name1)
        pose.dump_pdb(output_model_name)

# Protein Ligand concatenation


def protein_ligand_concatenation(template_hits, target_fasta_file, ligands):

    templates_df = pd.read_csv(template_hits, sep=",")
    # first templatePDB in a template Column
    top_first_hit_pdb = templates_df["template"].iloc[0]
    target_pdb_path = os.path.dirname(template_hits)
    target_pdb_name = os.path.basename(target_fasta_file).split('.')[0]
    top_hit_model = '{0}_{1}.pdb'.format(target_pdb_name, top_first_hit_pdb)

    for lig in ligands:
        protein_model = os.path.join(target_pdb_path, top_hit_model)
        basename_ligand = os.path.splitext(os.path.basename(lig))[0]
        complex_protein_ligand = '{0}_{1}.pdb'.format(
            top_hit_model.split('.')[0], basename_ligand)
        print('cat {0} {1} > protein_ligand_complex_top_1_comp_model/{2}'.format(
            protein_model, lig, complex_protein_ligand))
        os.system('cat {0} {1} > protein_ligand_complex_top_1_comp_model/{2}'.format(
            protein_model, lig, complex_protein_ligand))


currentWD = os.getcwd()
os.chdir(args['folder'])

# os.makedirs("protein_comp_modeling/protein_seq_alignment_files")
# emboss_needle_search(emboss_needle, glob.glob('*.fasta')[0], tempate_seqs)

# emboss_alignments = glob.glob(
#     'protein_comp_modeling/protein_seq_alignment_files/*.needle')
# rocs_single_rpt = 'ROCS/single_report_file_sorted.csv'
# select_top_hits_from_emboss_and_rocs_pdb(
#     emboss_alignments, rocs_single_rpt, glob.glob('*.fasta')[0])

template_hits = 'protein_comp_modeling/template_hits.csv'
# seq_alin_file = 'protein_comp_modeling/protein_seq_alignment_files/'
# modeling(template_hits, apo_pdbs, seq_alin_file, glob.glob('*.fasta')[0])

os.mkdir("protein_ligand_complex_top_1_comp_model")
protein_ligand_concatenation(template_hits, glob.glob(
    "*.fasta")[0], glob.glob("sdf2params/*.pdb"))

os.chdir(currentWD)
