function [min_energy] =  rotate_CH3_group_1_find_lowest(Position, HG1_Clash, protein_clash_HG1, rest_pro_position, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1)

done = 0;
Pos_b4_HG1 = Position;
min_loc = 0;
[HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
if HG1_energy == 0
    done = 1;
end
HG1_val = 0;
min_energy = HG1_energy;


                            
while done == 0 && HG1_val < 73
    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;
    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);

    if HG1_energy ==0
        done = 1;
        min_energy = HG1_energy;
        min_loc = setHG1;
    elseif HG1_energy < min_energy
        min_energy = HG1_energy;
        min_loc = setHG1;
    end

end
end
