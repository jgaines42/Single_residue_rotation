# Single_residue_rotation

Functions for rotating a single amino acid, either in a dipeptide or in the context of the protein. The main function is run_single_residue().

run_single_residue(PDB_name1, which_res,folder_name, save_folder, is_dipeptide)\
## Input:
PDB_name1: the name of the PDB to be run. Should be saved as XXXX_H.pdb or XXXX.mat where PDB_name1 = XXXX\
which_res: the residue ID\
folder_name: The full path name to the folder containing the PDB file.\
save_folder: Folder to save results to \
is_dipeptide: 1 to run just the dipeptide, 0 to run in context of the protein

## Notes:
Create PDB file with hydrogens by using download_preprocess_pdb.py which also adds the hydrogen atoms to the protein\
Final results are stored in \*protein_single_rotation_minE.mat or \*dipeptide_single_rotation_minE.mat

## Citation:
**reduce (used to add hydrogens):** Word, et al.(1999) "Asparagine and glutamine: using hydrogen atom contacts in the choice of sidechain amide orientation" J. Mol. Biol. 285, 1735-1747.\
**Single residue rotation:** J. C. Gaines, A. Virrueta, D. A. Buch, S. J. Fleishman, C. S. O'Hern, and L. Regan, "Collective repacking reveals that the structures of protein cores are uniquely specified by steric repulsive interactions," Protein Eng. Des. Sel. 30 (2017) 387.\
**Processing code:** J. C. Gaines, S. Acebes, A. Virrueta, M. Butler, L. Regan, and C. S. O'Hern, "Comparing side chain packing in soluble proteins, protein-protein interfaces, and transmembrane proteins," Proteins: Structure, Function, and Bioinformatics 86 (2018) 581. \
