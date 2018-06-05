%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [min_energy] =  rotate_CH3_group_2_find_lowest(Position, HG1_Clash, protein_clash_HG1, HG2_clash, protein_clash_HG2, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1,HG_Array_2, moveAtomID_HG2, numAtom,delta_term_HG2,global_min_sofar,E_sofar)
% 
% Finds the lowest energy configuration of 2 CH3 groups (or 1 CH3 and 1 OH group)
% resulting from rotating around the dihedral angles to the hydrogen atoms
%
% Input: 
% Position: Coordinates of the dipeptide
% HG1_Clash: Clash list within dipeptide to check with rotations of first
% group
% protein_clash_HG1: Clash list within protein for first group
% HG2_Clash: Clash list withing dipeptide to check with 2nd group
% protein_clash_HG2: Clash list with rest of protein for 2nd group
% rest_pro_postion: coordinates of the rest of the protein
% subtract_array_HG1: subtract array to put second atom of dihedral at
% orgin
% delta_term_HG1: delta term for Rotate_DA
% HG_Array_1: list of atoms defining the dihedral angle of group 1
% moveAtomID_HG1: atoms that move during rotation of group 1
% delta_term_HG2: delta term for Rotate_DA
% HG_Array_2: list of atoms defining the dihedral angle of group 2
% moveAtomID_HG2: atoms that move during rotation of group 2
%
% Output:
% min_energy: lowest energy of the combined rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [min_energy] =  rotate_CH3_group_2_find_lowest(Position, HG1_Clash, protein_clash_HG1, HG2_Clash, protein_clash_HG2, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1, HG_Array_2, moveAtomID_HG2, numAtom,delta_term_HG2,global_min_sofar,E_sofar)


Pos_b4_HG1 = Position;
HG1_val = 0;
done = 0;
min_overal = global_min_sofar+.01; %add 0.01 so that we can see if we found a configuration lower than the min_so_far

%%First check HG2 against the protein to make a short list of combinations to run
HG2_list = [];
    ind0 = ismember(HG2_Clash(:,1), [1:8]);

for i = 1:72
    setHG2 = i*5;
    Position=Pos_b4_HG1;
    subtract_array_HG2 = repmat(Position(HG_Array_2(2),:),numAtom,1);
    Position = Rotate_DA(Position, setHG2, subtract_array_HG2, delta_term_HG2, HG_Array_2, moveAtomID_HG2);
    if size(protein_clash_HG2,1) > 0
        [HG2_energy] = get_energy_wProtein([], Position,  protein_clash_HG2, rest_pro_position);
    else
        HG2_energy =0;
    end
    % Check with part of backbone that is the same with every amino acid
    HG2_energy = HG2_energy + get_energy_wProtein(HG2_Clash(ind0,:), Position,[], rest_pro_position);

    if HG2_energy + E_sofar <= global_min_sofar
        HG2_list = [HG2_list,i];
    end
end

%If none of HG2 rotations had low energy with the rest of the protein,
%return to prior function
if size(HG2_list,1) == 0
    done = 1;
    min_energy = global_min_sofar + 1;
end

%loop over all HG1 values and sample HG2 till low energy state is found
while done == 0 && HG1_val <= 72

    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position= Pos_b4_HG1;

    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
 
    if HG1_energy >0
        no_Clash = 0; %there are clashes
    else
        no_Clash = 1; %no HG1 clashes
    end
    
    if (HG1_energy < min_overal && HG1_energy + E_sofar < global_min_sofar) || no_Clash == 1 %if HG1 energy is low or no clashes so far
        Pos_b4_HG2 = Position;
        HG2_val = 0;
        subtract_array_HG2 = repmat(Position(HG_Array_2(2),:),numAtom,1);
        no_2_clash = 0;
        
        while done == 0 && HG2_val < size(HG2_list,2) && no_2_clash == 0
            
            HG2_val = HG2_val +1;
            setHG2 = HG2_list(HG2_val)*5;
            Position = Pos_b4_HG2;
            
            Position = Rotate_DA(Position, setHG2, subtract_array_HG2, delta_term_HG2, HG_Array_2, moveAtomID_HG2);
            [HG2_energy] = get_energy_wProtein(HG2_Clash, Position,  protein_clash_HG2, rest_pro_position);
            if HG2_energy >0 && HG2_val == 1
                min_HG2 = HG2_energy;
                no_2_clash = 0;
            elseif HG2_energy >0 && HG2_energy < min_HG2
                min_HG2 = HG2_energy;
                no_2_clash = 0;
            elseif HG2_energy == 0
                no_2_clash = 1;
                min_HG2 = 0;
                if no_Clash == 1
                    done = 1;
                end
            end
        end
        
        if (HG1_energy + min_HG2) < min_overal
            min_overal = HG1_energy + min_HG2;

          
        end
    end
end
min_energy = min_overal;
end
