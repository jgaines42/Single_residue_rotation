
% total_energy must be the energy so far for the residue
% energy_cuttoff is set in the main function. All CH3 movement that result
% in a total energy less than this will be kept

done = 0;
degrees_to_sample = 72/3; %cause there is rotational symmetry

Pos_b4_HG1 = Position;

HG1_val = 0;
while  HG1_val < degrees_to_sample %done == 0 &&
    HG1_val = HG1_val +1;
    setHG1 = HG1_val*5;
    Position=Pos_b4_HG1;
    Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
    [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);

    if HG1_energy + total_energy <= energy_cutoff        
        energy(count_lowE,1) = total_energy + HG1_energy;
        coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
        dihedrals(count_lowE,1:4) =[setChi1, setChi2, setChi3, setHG1];
        count_lowE = count_lowE+1;
        
        if size(dihedrals,1) == (count_lowE -1)
            energy = [energy; zeros(10000,1)];
            coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
            dihedrals = [dihedrals;zeros(10000,DOF+CH3)];
        end
                
    end

end
