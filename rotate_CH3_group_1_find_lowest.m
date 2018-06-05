%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [min_energy] =  rotate_CH3_group_1_find_lowest(Position, HG1_Clash, protein_clash_HG1, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1)
% 
% Finds the lowest energy configuration of a CH3 group (or OH group)
% resulting from rotating around the dihedral angle to the hydrogen atom
%
% Input: 
% Position: Coordinates of the dipeptide
% HG1_Clash: Clash list within dipeptide to check with rotations
% protein_clash_HG1: Clash list within protein
% rest_pro_postion: coordinates of the rest of the protein
% subtract_array_HG1: subtract array to put second atom of dihedral at
% orgin
% delta_term_HG1: delta term for Rotate_DA
% HG_Array_1: list of atoms defining the dihedral angle
% moveAtomID_HG1: atoms that move during rotation
%
% Output:
% min_energy: lowest energy of the rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [min_energy] =  rotate_CH3_group_1_find_lowest(Position, HG1_Clash, protein_clash_HG1, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1)


Pos_b4_HG1 = Position; % Store original position

%Check to see if original position has energy of 0
[HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
if HG1_energy == 0
    done = 1;
end
HG1_val = 0;
min_energy = HG1_energy;
done = 0;



% While we haven't found a spot with energy== 0, loop through angles                            
while done == 0 && HG1_val < 73
    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;
    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);

    if HG1_energy ==0
        done = 1;
        min_energy = HG1_energy;
    elseif HG1_energy < min_energy
        min_energy = HG1_energy;
    end

end
end
