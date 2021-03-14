# Rosetta Kinase CM

<!-- <figure class="image">
  <img src="pipeline.png">
  <figcaption>Figure: Overview of Rosetta Kinase CM pipeline</figcaption>
</figure> -->

## Dependency
* Tested with [Python 3.7](https://www.python.org/downloads/) (Pandas) 
* [Rosetta](https://www.rosettacommons.org/software/license-and-download) Software Suite 
* [PyRosetta](http://www.pyrosetta.org/) Software Suite
* [OpenEye](https://www.eyesopen.com/) Software Suite
* [EMBOSS](http://emboss.open-bio.org/html/use/ch02s07.html) Software Suite

Initial folder should contain these three files and similar naming convention,

```
2W1C_A_L0C
├── 2W1C_A.fasta
├── 2W1C_A_L0C.pdb
└── L0C.smi
```
Here are the STEPS and SCRIPTS used for the comparative modeling pipeline approach,

## 1. Conformer generation and ligand alignment (modeling_script.py)
* First, the maximum number of conformers will be generated for a given SMILES using OMEGA. The output will be a single SDF file with the maximum number of conformers.
* Second, the single SDF file will be used for the alignment of a query (Molecule of interest) and database (in-house active kinase ligand template library) molecules using ROCS. The output will be each template aligned query molecule.
* Third, combining the report files of each template aligned query molecule into a single report file.
* Fourth, based on the single report file, the top 100 conformers for the given SMILES are selected.
* Fifth, SDF to PARAMS will be generated using 100 conformers for Rosetta Minimization step. 
```
python /path/to/Rosetta_Kinase_CM/modeling_script.py -f /path/to/2W1C_A_L0C -omega /path/to/omega2 -rocs /path/to/rocs -temp_lig /path/to/Rosetta_Kinase_CM/template_ligand_library -mol2params /path/to/Rosetta/main/source/scripts/python/public/generic_potential/mol2genparams.py -convert /path/to/Rosetta_Kinase_CM/convert.py
```

## Output folder should contain similar files after running the above command

```
2W1C_A_L0C
├── 2W1C_A.fasta
├── 2W1C_A_L0C.pdb
├── L0C.smi
├── OMEGA
│   ├── L0C_omega.log
│   ├── L0C_omega.parm
│   ├── L0C_omega.rpt
│   ├── L0C_omega.sdf
│   └── L0C_omega_status.txt
├── ROCS [7787 entries exceeds filelimit, not opening dir]
├── mol2params [200 entries exceeds filelimit, not opening dir]
└── top_100_conf [200 entries exceeds filelimit, not opening dir]
```

## 2. Sequence alignment and protein modeling (new_protein_modeling.py)
* First, target-template sequence alignment will be performed using in-house active kinase sequence template library (EMBOSS Needleman-Wunsch algorithm)
* Second, selection of top hit templates will be applied using sequence and ligand similarity approach (defined as Template Score). Based on this approach, the top 10 templates for the given sequence will be selected.
* Third, 10 predicted models of target protein will be performed using PyRosetta.
* Forth, using the top first model we concatenate the 100 conformers from ligand alignment and this results into an unrefined protein-ligand complex of 100 comparative models.
```
python new_protein_modeling.py -f /path/to/2W1C_A_L0C -emboss /path/to/usr/local/emboss/bin/needle -temp_seq /path/to/rosetta_kinase_cm/template_fasta_seq_training_set -apo_pdb /path/to/rosetta_kinase_cm/apo_pdbs_for_template_seq_extraction
```

## Output folder should contain similar files after running the above command
```
2W1C_A_L0C
├── 2W1C_A.fasta
├── 2W1C_A_L0C.pdb
├── L0C.smi
├── OMEGA
│   ├── L0C_omega.log
│   ├── L0C_omega.parm
│   ├── L0C_omega.rpt
│   ├── L0C_omega.sdf
│   └── L0C_omega_status.txt
├── ROCS [7787 entries exceeds filelimit, not opening dir]
├── mol2params [200 entries exceeds filelimit, not opening dir]
├── protein_comp_modeling [13 entries exceeds filelimit, not opening dir]
├── protein_ligand_complex_top_1_comp_model [100 entries exceeds filelimit, not opening dir]
└── top_100_conf [200 entries exceeds filelimit, not opening dir]
```

## 3. Minimization of protein-ligand complex (minimization.py)
* First, the input files for Rosetta minimization process will be generated for parallel computing.
* Second, once minimization finished, the energy for each model will be calculated similarly to the first step.
## 4. Analysis (analysis_1.py)
* Here, I will generate the table for 100 minimized structures. It contains the name and energy attributes of those models. Out of 100, the top 10 models will be selected using Rosetta energy values.
## 5. Complex modeling of remaining protein models (top_comp_prtn_lig_modeling.py) (from step 2, third point)
* Here, the PARAMS files of top 10 ligands will be taken from step 1, fifth point.
* Concatenation of protein-ligand complex (this will again result into 100 complex models) 
## 6. Minimization (top_comp_prtn_lig_modeling_minimization.py)
* The minimization process is same as step 3.
## 7. Analysis (analysis_2.py)
* The analysis process is same as step 4.
* The top 1 model will be reported as the best prediction. 

# Contact
Reach me at kirubanpk@gmail.com

# License
This project uses the following license: MIT License
