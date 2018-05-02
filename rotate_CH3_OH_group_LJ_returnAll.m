Pos_b4_HG1 = Position;
subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;


HG1_val = 0;
degrees_to_sample = 72/3; %cause there is rotational symmetry

while  HG1_val < degrees_to_sample
    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;
    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
    
    
    if HG1_energy +total_energy < energy_cutoff
        Pos_b4_OH = Position;
        subtract_array_OH = repmat(Position(iOHArray(2),:),numAtom,1);
        delta_term_OH =  pi*sign(InitOH)*InitOH/180;
        
        OH_val = 0;
        while OH_val < 72
            OH_val = OH_val +1;
            setOH = OH_val*5;
            Position = Pos_b4_OH;
            Position = Rotate_DA(Position, setOH, subtract_array_OH, delta_term_OH, iOHArray, moveAtomOH);
            [OH_energy] = get_energy_wProtein(OH_Clash, Position,  protein_clash_OH, rest_pro_position);
            
            if OH_energy + HG1_energy + total_energy < energy_cutoff
                
                energy(count_lowE,1) = total_energy + OH_energy + HG1_energy;
                coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                dihedrals(count_lowE,1:3) = [setChi1,setHG1,setOH];
                count_lowE = count_lowE+1;
                if size(dihedrals,1) == (count_lowE -1)
                    energy = [energy; zeros(10000,1)];
                    coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                    dihedrals = [dihedrals;zeros(10000,DOF+CH3+OH)];
                end
            end
                
        end
    end
end