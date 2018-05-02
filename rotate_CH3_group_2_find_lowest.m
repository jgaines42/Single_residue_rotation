function [min_energy] =  rotate_CH3_group_2_find_lowest(Position, HG1_Clash, protein_clash_HG1, HG2_Clash, protein_clash_HG2, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1, HG_Array_2, moveAtomID_HG2, numAtom,delta_term_HG2,global_min_sofar,E_sofar)


Pos_b4_HG1 = Position;
HG1_val = 0;
done = 0;
min_overal = global_min_sofar; %min energy
gmin_loc = [0,0];

%%First check HG2 against the protein to make a short list of combinations to run
HG2_list = [];
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
    if HG2_energy + E_sofar <= global_min_sofar
        HG2_list = [HG2_list,i];
    end
end
if size(HG2_list,1) == 0
    done = 1;
    min_energy = global_min_sofar + 1;
end

%loop over all HG1 values
while done == 0 && HG1_val <= 72

    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;

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
