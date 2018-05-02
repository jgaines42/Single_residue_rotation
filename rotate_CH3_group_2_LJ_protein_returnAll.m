% total_energy must be the energy so far for the residue
% energy_cuttoff is set in the main function. All CH3 movement that result
% in a total energy less than this will be kept

degrees_to_sample = 72/3; %cause there is rotational symmetry
Pos_b4_HG1 = Position;
HG1_val = 0;
done = 0;

%loop over all HG1 values
while  HG1_val < degrees_to_sample

    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;
    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
  
    if (HG1_energy+total_energy <= energy_cutoff) 
        Pos_b4_HG2 = Position;
        HG2_val = 0;
        subtract_array_HG2 = repmat(Position(HG_Array_2(2),:),numAtom,1);
        delta_term_HG2 =  pi*sign(InitChi_HG2)*InitChi_HG2/180;

        while  HG2_val < degrees_to_sample 
            
            HG2_val = HG2_val +1;
            setHG2 = HG2_val*5;
            Position = Pos_b4_HG2;
            
            Position = Rotate_DA(Position, setHG2, subtract_array_HG2, delta_term_HG2, HG_Array_2, moveAtomID_HG2);
            [HG2_energy] = get_energy_wProtein(HG2_Clash, Position,  protein_clash_HG2, rest_pro_position);
           
            if (HG1_energy + HG2_energy + total_energy <= energy_cutoff)
                energy(count_lowE,1) = total_energy + HG1_energy + HG2_energy;
                coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                if DOF == 1
                     dihedrals(count_lowE,1:3) =[setChi1, setHG1, setHG2];
                elseif DOF == 2
                     dihedrals(count_lowE,1:4) =[setChi1, setChi2, setHG1, setHG2];
                else
                    dihedrals(count_lowE,1:5) =[setChi1, setChi2, setChi3, setHG1, setHG2];
                end
                count_lowE = count_lowE+1;
                if size(dihedrals,1) == (count_lowE -1)
                    energy = [energy; zeros(10000,1)];
                    coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                    dihedrals = [dihedrals;zeros(10000,DOF+CH3)];
                end
            end

        end
    end
end