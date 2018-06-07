%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [total_energy] = get_energy_of_dipeptide(Clash_lists, Position, DOF, CH3, OH)
%
% Calculates the total energy of a given dipeptide
%
% Input:
% Clash_lists: Structure containing all the clash lists for the amino acid
% Position: Positions of the atoms
% DOF: Number of regular side chain dihedral angles
% CH3: Number of CH3 groups
% OH: Number of OH groups
%
% Returns:
% total_energy: The total repulsive Lennard-Jones energy of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [total_energy] = get_original_energy(Clash_lists, Position, DOF, CH3, OH,rest_of_pro)

total_energy = 0;
if DOF >= 1
    
    [chi1_energy] = get_energy_wProtein(Clash_lists.c1_Clash, Position, Clash_lists.protein_clash_c1, rest_of_pro);
    total_energy = total_energy +chi1_energy;
    
    if DOF >= 2
        [chi2_energy] = get_energy_wProtein(Clash_lists.c2_Clash, Position,   Clash_lists.protein_clash_c2, rest_of_pro);
        total_energy = total_energy + chi2_energy;
    end
    if DOF >= 3
        [chi3_energy] = get_energy_wProtein(Clash_lists.c3_Clash, Position,  Clash_lists.protein_clash_c3, rest_of_pro);
        total_energy = total_energy + chi3_energy;
    end
end
if CH3 == 2
    [HG2_energy] = get_energy_wProtein(Clash_lists.HG2_Clash, Position,  Clash_lists.protein_clash_HG2, rest_of_pro);
    total_energy = total_energy + HG2_energy;
end
if CH3 >=1
    [HG1_energy] = get_energy_wProtein(Clash_lists.HG1_Clash, Position,  Clash_lists.protein_clash_HG1, rest_of_pro);
    total_energy = total_energy + HG1_energy;
end
if OH ==1
    [OH_energy] = get_energy_wProtein(Clash_lists.OH_Clash, Position,   Clash_lists.protein_clash_OH, rest_of_pro);
    total_energy = total_energy + OH_energy;
end
end