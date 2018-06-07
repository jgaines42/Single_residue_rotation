%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[dihedrals, min_energy] = check_end_groups(dihedrals, min_energy,total_energy, this_dihedral, Amino_Acid,Init_dihedrals, Position, Clash_list)
%
% Determines the lowest energy of the system after rotating the OH or CH3
% groups at the end of the protein
%
% Input:
% dihedrals: List of dihedrals so far that have min_energy
% min_energy: global minimum energy so far for the amino acid
% this_diheral: chi1, chi2, etc for this configuration
% Amino_acid: Structure containing rotation data for the amino acid type
% Init_dihedrals: Structuer containing initial dihedral angles for the amino aicd
% Position: coordinates of the amino acid
% Clash_list: Structure containing clash lists for the amino acid
%
% Output:
% dihedrals: Array of all dihedral angles with min_energy
% min_energy: global minimum energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dihedrals, min_energy] = check_end_groups(dihedrals, min_energy,total_energy, this_dihedral, Amino_Acid,Init_dihedrals, Position, Clash_list,rest_of_pro)

% Extract data from Structures
CH3 = Amino_Acid.CH3;
OH = Amino_Acid.OH;
numAtom = Amino_Acid.numAtom;


if CH3 == 1 &&  OH == 0 % Met
    subtract_array_HG1 = repmat(Position(Amino_Acid.HG_Array_1(2),:),numAtom,1);
    delta_term_HG1 =  pi*sign(Init_dihedrals.InitChi_HG1)*Init_dihedrals.InitChi_HG1/180;
    min_energy_CH3  = rotate_CH3_group_1_find_lowest(Position, Clash_list.HG1_Clash, Clash_list.protein_clash_HG1, rest_of_pro, subtract_array_HG1, delta_term_HG1, Amino_Acid.HG_Array_1, Amino_Acid.moveAtomID_HG1);   
    total_energy = total_energy + min_energy_CH3;

elseif CH3 == 2 &&  OH == 0 % Val, Ile, Leu, 
    
    subtract_array_HG1 = repmat(Position(Amino_Acid.HG_Array_1(2),:),numAtom,1);
    delta_term_HG1 =  pi*sign(Init_dihedrals.InitChi_HG1)*Init_dihedrals.InitChi_HG1/180;
    delta_term_HG2 =  pi*sign(Init_dihedrals.InitChi_HG2)*Init_dihedrals.InitChi_HG2/180;
    min_energy_CH3 = rotate_CH3_group_2_find_lowest(Position, Clash_list.HG1_Clash, Clash_list.protein_clash_HG1, Clash_list.HG2_Clash, Clash_list.protein_clash_HG2, rest_of_pro, subtract_array_HG1, delta_term_HG1, Amino_Acid.HG_Array_1, Amino_Acid.moveAtomID_HG1, Amino_Acid.HG_Array_2, Amino_Acid.moveAtomID_HG2, numAtom,delta_term_HG2,min_energy, total_energy);
    total_energy = total_energy + min_energy_CH3;

elseif CH3 == 0 && OH == 1 %Ser, Tyr
    Pos_b4_OH = Position;
    subtract_array_OH = repmat(Position(Amino_Acid.iOHArray(2),:),numAtom,1);
    delta_term_OH =  pi*sign(Init_dihedrals.InitOH)*Init_dihedrals.InitOH/180;
    min_energy_OH  = rotate_CH3_group_1_find_lowest(Position, Clash_list.OH_Clash, Clash_list.protein_clash_OH, rest_of_pro, subtract_array_OH, delta_term_OH, Amino_Acid.iOHArray, Amino_Acid.moveAtomOH);
    total_energy = total_energy + min_energy_OH;
    
elseif CH3 == 1 && OH == 1 % Thr
    Pos_b4_HG1 = Position;
    
    subtract_array_HG1 = repmat(Position(Amino_Acid.HG_Array_1(2),:),numAtom,1);
    delta_term_HG1 =  pi*sign(Init_dihedrals.InitChi_HG1)*Init_dihedrals.InitChi_HG1/180;   
    delta_term_OH = pi*sign(Init_dihedrals.InitOH)*Init_dihedrals.InitOH/180;
    min_energy_1 = rotate_CH3_group_2_find_lowest(Position, Clash_list.HG1_Clash, Clash_list.protein_clash_HG1, Clash_list.OH_Clash, Clash_list.protein_clash_OH, rest_of_pro, subtract_array_HG1, delta_term_HG1, Amino_Acid.HG_Array_1, Amino_Acid.moveAtomID_HG1, Amino_Acid.iOHArray, Amino_Acid.moveAtomOH, numAtom,delta_term_OH,min_energy, total_energy);
    total_energy = total_energy + min_energy_1;
    
    
end

% If total_energy of this configuration of the amino acid is lower than
% anything found so far, clear dihedrals and just return this amino acids
% configuration and change min_energy
if total_energy < min_energy
    min_energy = total_energy;
    dihedrals = this_dihedral;
elseif total_energy == min_energy %If total_energy is the same, just concatinate the diheral
    dihedrals = [dihedrals;this_dihedral];
end



end