
function [using_E] = find_minimum_energy_toC_bash_noCH3(pdb_names, variant,database, group, min_tostart)
tic
 if strcmp(database,'i')
    folder1 = strcat('Interface_rotations/', pdb_names, '/');
    folder_name = '/home2/jcg72/volumes/Interface_database/';
elseif strcmp(database,'h')
  % folder1 = strcat('HQ_54_multi/', pdb_names, '/');
  %  folder_name = '/Users/jennifergaines/Documents/summer2013/Minimum_energy_algorithm/Core_mutations/';
  %  folder_name = '/home2/jcg72/Minimum_energy/HQ_input/';
  %  folder1 = strcat('/fastscratch/jcg72/HQ_54_multi/', pdb_names(1:4),'/');
    folder_name = '/home/fas/ohern/jcg72/Minimum_energy/HQ_input/';
    folder1 = strcat('/scratch/fas/ohern/jcg72/HQ_54_multi/',pdb_names, '/');
elseif strcmp(database,'m')
    folder1 = strcat('Core_mutations/',pdb_names(1:4), '/');
    folder_name = '/Users/jennifergaines/Documents/summer2013/Minimum_energy_algorithm/Core_mutations/';
    %folder_name = '/home2/jcg72/Minimum_energy/Core_mutations/';
elseif strcmp(database,'s')
    folder1 = '/scratch/fas/ohern/jcg72/STSC_IW/';
    %folder1 = '/home2/jcg72/STSC_IW/';
    %folder_name = folder1;
    %folder1 = '/Users/jennifergaines/Documents/IW_Fall2016/';
    %folder_name = '/Users/jennifergaines/Documents/IW_Fall2016/';

end


if min_tostart > 1
    f = fopen(strcat(folder1, pdb_names, '_loc_out_', num2str(variant),'_g', num2str(group), '.txt'));
    using_E = fscanf(f,'%f');
    fclose(f);
    using_E = using_E*2;
    min_tostart = using_E;
else
    using_E = 0.0001;
    min_tostart = 0.0001;
end
 
if ~exist(strcat(folder1,pdb_names, '_parameters_', num2str(variant),'_g', num2str(group),  '_doAll.txt'), 'file')
    
    if strcmp(database,'i')
        PDB_name = strcat(folder_name, pdb_names, '_ordered.mat');
        load(PDB_name);
        load(strcat(folder_name, 'Interface_9sizes/', pdb_names, '_interface_residues_new.mat'));
        all_moving = double(inter_core);
        
    elseif strcmp(database, 's')
        PDB_name = strcat(folder_name, pdb_names, '.mat');
        load(PDB_name);
        load(strcat(folder_name, 'STSC_moving.mat'));
        ind0 = new_moving(:,2) == group;
        all_moving = double(new_moving(ind0,1));
    else
        PDB_name = strcat(folder_name, pdb_names, '.mat');
        load(PDB_name);
        load(strcat(folder_name, pdb_names, '_all_core_binned_2.mat'));
        ind0 = new_moving(:,2) == group;
        all_moving = double(new_moving(ind0,1));
    end
    
    ind1 = strcmp(tempModel2(:,3),'B');
    x = tempModel2(~ind1,:);
    tempModel2 = tempModel2(~ind1,:);
    tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
    
    %relabel residue ids of second chain so all > ids of 1st chain
    res_ids = cell2mat(tempModel2(:,6));
    for i = 2:size(tempModel2,1)
        if res_ids(i,1) < res_ids(i-1,1)
            res_ids(i,1) = res_ids(i,1) + 1000;
        end
    end
    tempModel2(:,6) = num2cell(res_ids);
    tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
    
    %Add sizes
    tempModel2 = add_sizes_protein(tempModel2,9);
    sizes_all = cell2mat(tempModel2(:,12));
    
    Residue_id = [];
    Residue_names = {};
    all_DOF = zeros(size(all_moving,1),1);
    for i = 1:size(all_moving,1)
        ind0 = find(res_ids == all_moving(i));
        resiName = tempModel2{ind0(1),4}; %Get name of interface resiude
        resiName(2:3) = lower(resiName(2:3));
        if database == 'h'
            Residue_id = [Residue_id,(all_moving(i))];
            Residue_names = [Residue_names, resiName];
        elseif ismember(resiName, { 'Ala','Ile', 'Leu', 'Phe', 'Met', 'Val'})
            Residue_id = [Residue_id,(all_moving(i))];
            Residue_names = [Residue_names, resiName];
        end
     
        
        if  resiName == 'Leu' | resiName == 'Ile' | resiName == 'Phe' | resiName == 'His' | resiName == 'Tyr' | resiName == 'Trp'
            all_DOF(i) = 2;
        elseif resiName == 'Met'
            all_DOF(i) = 3;
        else
            all_DOF(i) = 1;
        end
    end
    all_moving = Residue_id;
    
    
    % Get the right varients of each
    for i = 1:size(all_moving,2)
        % Isolate the dipeptide and make sure it has the right order of atoms
        [allDipeptide,next_pro] = isolate_dipeptide(tempModel2, res_ids, all_moving(i));
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-3,:), Residue_names{i});
        allDipeptide(4:size(allDipeptide,1)-3,:)= new_Dipeptide;
        ind0 = ismember(cell2mat(tempModel2(:,1)), cell2mat(allDipeptide(:,1)));
        tempModel2(ind0,:) = allDipeptide;
        tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
        if next_pro
            allDipeptide = [allDipeptide;num2cell(zeros(1, size(allDipeptide,2)))];
        end
        
        ind0 = ismember(cell2mat(tempModel2(:,1)), cell2mat(allDipeptide(:,1)));
        rest_of_pro = tempModel2(~ind0,:);
        
        [Position, rest_of_pro] =get_variant(allDipeptide, Residue_names{i}, rest_of_pro, variant);
        if next_pro
            tempModel2(ind0,8:10) = num2cell(Position(1:size(Position,1)-1,:));
        else
            tempModel2(ind0,8:10) = num2cell(Position);
        end
        tempModel2(~ind0,:) = rest_of_pro;
        
    end
    
    %% Now the the protein is set up, find the lowest energy
    
    if size(Residue_id,2)>1 && size(Residue_id,2) < 11
        
        
        
        min_found = 0; %Have we found the minimum yet
        
        clear coord_1;
        clear coord_2;
        clear coord_3;
        clear coord_4;
        clear coord_5;
        clear coord_6;
        clear coord_7;
        clear coord_8;
        clear coord_9;
        clear coord_10;
        
        global coord_1 coord_2 coord_3 coord_4 coord_5 coord_6 coord_7 coord_8 coord_9 coord_10;
        coor_m1 = [];
        coor_m2 = [];
        coor_m3 = [];
        coor_m4 = [];
        coor_m5 = [];
        coor_m6 = [];
        coor_m7 = [];
        coor_m8 = [];
        coor_m9 = [];
        coor_m10 = [];
        
        residues = 1; %Loop over all residues
        continue_res = 1; %Should we test the next residue?
        mins = zeros(1,size(Residue_id,2)); %All local min E values
        not_found = 1;
        
        while not_found
            residues = 1;
            continue_res = 1;
            while residues <= size(Residue_id,2) && continue_res;
                
                [all_coordinates, dihedrals] = set_up_recursion_SingleRes_rot_noCH3(Residue_id(residues), Residue_names{residues}, tempModel2, all_moving, min_tostart);
                if size(dihedrals,1) > 0
                    save(strcat(folder1, pdb_names, '_', Residue_names{residues}, num2str(Residue_id(residues)), '_dihedrals_fast2_', num2str(variant), '.mat'), 'dihedrals');
                    switch residues
                        case 1
                            coord_1 = all_coordinates;
                            energy1 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            [i,loc] = min(energy1);
                            mins(residues) = i;
                            coor_m1 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 2
                            coord_2 = all_coordinates;
                            energy2 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy2);
                            [i,loc] = min(energy2);
                            mins(residues) = i;
                            coor_m2 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                        case 3
                            coord_3 = all_coordinates;
                            energy3 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy3);
                            [i,loc] = min(energy3);
                            mins(residues) = i;
                            coor_m3 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 4
                            coord_4 = all_coordinates;
                            energy4 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy4);
                            
                            [i,loc] = min(energy4);
                            mins(residues) = i;
                            coor_m4 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                        case 5
                            coord_5 = all_coordinates;
                            energy5 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy5);
                            [i,loc] = min(energy5);
                            mins(residues) = i;
                            coor_m5 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 6
                            coord_6 = all_coordinates;
                            energy6 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy6);
                            [i,loc] = min(energy6);
                            mins(residues) = i;
                            coor_m6 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 7
                            coord_7 = all_coordinates;
                            energy7 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy7);
                            [i,loc] = min(energy7);
                            mins(residues) = i;
                            coor_m7 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 8
                            coord_8 = all_coordinates;
                            energy8 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy8);
                            [i,loc] = min(energy8);
                            mins(residues) = i;
                            coor_m8 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 9
                            coord_9 = all_coordinates;
                            energy9 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy9);
                            [i,loc] = min(energy9);
                            mins(residues) = i;
                            coor_m9 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                        case 10
                            coord_10 = all_coordinates;
                            energy10 = all_coordinates(size(all_coordinates,1), 1:3:size(all_coordinates,2)-1);
                            mins(residues) = min(energy10);
                            [i,loc] = min(energy10);
                            mins(residues) = i;
                            coor_m10 = all_coordinates(:,[(loc-1)*3+1:loc*3,size(all_coordinates,2)]);
                            
                    end
                else % If nothing was returned, exit the loop
                    continue_res = 0;
                end
                
                residues = residues + 1;
            end
            if continue_res == 1
                not_found = 0;
            else
                if min_tostart < 0.01
                 min_tostart = min_tostart*5;
                else
                    min_tostart = min_tostart*2;
                end
            end
        
        
        
        using_E = min_tostart;
        %After checking all of the residues, start finding global min
        if continue_res
            min_loc = zeros(1,size(Residue_id,2));
            min_energy = min_tostart;
            max_depth = size(Residue_id,2);
            
            %See what the energy is of the combo of the lowest energy states
            if max_depth <= 10
                min_energy_1 = set_energy_cutoff_2(mins,coor_m1,coor_m2,coor_m3,coor_m4,coor_m5,coor_m6,coor_m7, coor_m8, coor_m9, coor_m10,max_depth);
                if min_energy_1 < min_energy
                    min_energy = min_energy_1;
                    for resi = 1:size(Residue_id,2)
                        current_coord = eval(strcat('coord_', num2str(resi)));
                        current_e =eval(strcat('energy', num2str(resi)));
                        load(strcat(folder1, pdb_names, '_', Residue_names{resi}, num2str(Residue_id(resi)), '_dihedrals_fast2_', num2str(variant), '.mat'));
                        ind1 = find(current_e <= min_energy);
                        ind2 = sort([(ind1-1)*3+1, (ind1-1)*3+2, (ind1-1)*3+3]);
                        dihedrals = dihedrals(ind1,:);
                        if size(dihedrals,1) == 0
                            not_found = 1;
                        end
                        save(strcat(folder1, pdb_names, '_', Residue_names{resi}, num2str(Residue_id(resi)), '_dihedrals_fast2_', num2str(variant), '.mat'), 'dihedrals');
                        
                        switch resi
                            case 1
                                coord_1 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy1 = current_e(ind1);
                            case 2
                                coord_2 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy2 = current_e(ind1);
                            case 3
                                coord_3 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy3 = current_e(ind1);
                            case 4
                                coord_4 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy4 = current_e(ind1);
                            case 5
                                coord_5 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy5 = current_e(ind1);
                            case 6
                                coord_6 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy6 = current_e(ind1);
                            case 7
                                coord_7 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy7 = current_e(ind1);
                            case 8
                                coord_8 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy8 = current_e(ind1);
                            case 9
                                coord_9 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy9 = current_e(ind1);
                            case 10
                                coord_10 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                energy10 = current_e(ind1);
                        end
                        
                    end
                end
                
            end
            
            if max_depth > 2
                %Do pairwise of lowest depth
                min_sofar = min_energy;
                
                %Do random sampling
                [overall_min,best_loc] = set_energy_cutoff_3( max_depth, min_energy);
                
                if overall_min < min_energy %If a low global minimum was found, change all the data to match this
                    min_energy = overall_min;
                    min_sofar = min_energy;
                    
                    for resi = 1:size(Residue_id,2)
                        current_coord = eval(strcat('coord_', num2str(resi)));
                        current_e =eval(strcat('energy', num2str(resi)));
                        load(strcat(folder1, pdb_names, '_', Residue_names{resi}, num2str(Residue_id(resi)), '_dihedrals_fast2_', num2str(variant), '.mat'));
                        ind1 = find(current_e <= min_energy);
                        ind2 = sort([(ind1-1)*3+1, (ind1-1)*3+2, (ind1-1)*3+3]);
                        dihedrals = dihedrals(ind1,:);
                        save(strcat(folder1, pdb_names, '_', Residue_names{resi}, num2str(Residue_id(resi)), '_dihedrals_fast2_', num2str(variant), '.mat'), 'dihedrals');
                        
                        switch resi
                            case 1
                                coord_1 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 2
                                coord_2 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                                
                            case 3
                                coord_3 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 4
                                coord_4 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 5
                                coord_5 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 6
                                coord_6 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 7
                                coord_7 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 8
                                coord_8 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 9
                                coord_9 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                            case 10
                                coord_10 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                        end
                        
                    end
                end
                
                
                
               %%%%%%%%%%%%%%%%%%
               % Run all pairwise based on the connectivity plot
                % run_pairwise;
               %%%%%%%%%%%%%%%%%%
                  here = 1; 
                
            end
            %  tic
            % min_loc = zeros(1000,max_depth);
            %[min_sofar, loc_min] = E_recursive_mins(0, min_energy, min_loc,[], 1, max_depth, zeros(1,size(Residue_id,2)),mins);
            % if min_energy < min_tostart %Write to file so you can run c++ code
            %max(max(loc_min)) > 0
           
            if not_found == 0
            for to_write = 1:max_depth
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_Xcoordinates_', num2str(variant), '.txt'), 'w');
                current_coord = eval(strcat('coord_', num2str(to_write)));
                for dihedr = 1:(size(current_coord,2)-1)/3
                    for atoms = 1:size(current_coord,1)-1
                        fprintf(f, '%f \n', current_coord(atoms,(dihedr-1)*3+1));
                    end
                end
                fclose(f);
                
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_Ycoordinates_', num2str(variant), '.txt'), 'w');
                current_coord = eval(strcat('coord_', num2str(to_write)));
                for dihedr = 1:(size(current_coord,2)-1)/3
                    for atoms = 1:size(current_coord,1)-1
                        fprintf(f, '%f \n', current_coord(atoms,(dihedr-1)*3+2));
                    end
                end
                fclose(f);
                
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_Zcoordinates_', num2str(variant), '.txt'), 'w');
                current_coord = eval(strcat('coord_', num2str(to_write)));
                for dihedr = 1:(size(current_coord,2)-1)/3
                    for atoms = 1:size(current_coord,1)-1
                        fprintf(f, '%f \n', current_coord(atoms,(dihedr-1)*3+3));
                    end
                end
                fclose(f);
                
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_energy_', num2str(variant), '.txt'), 'w');
                energyX = current_coord(size(current_coord,1), [1:3:size(current_coord,2)-1]);
                for dihedrals = 1:size(energyX,2)
                    fprintf(f, '%f \n', energyX(dihedrals));
                end
                fclose(f);
                
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_atomSizes_', num2str(variant), '.txt'), 'w');
                atomSizes = current_coord(1:size(current_coord,1)-1, size(current_coord,2));
                for atoms = 1:size(atomSizes,1)
                    fprintf(f, '%f \n', atomSizes(atoms));
                end
                fclose(f);
                
               %Assign everything to unique groups (cause no CH3 stuff
               load(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_dihedrals_fast2_', num2str(variant), '.mat'));
                this_DOF = all_DOF(to_write);
                if this_DOF < size(dihedrals,2)
                    unique_dih = unique(dihedrals(:,1:this_DOF), 'rows');
                    
                    dih_groups = zeros(size(dihedrals,1),1);
                    start_dih = dihedrals(1,1:this_DOF);
                    dih_groups(1) = 1;
                    group_num = 1;
                    for d = 2:size(dihedrals,1)
                        if (dihedrals(d,1:this_DOF) == start_dih)
                            dih_groups(d) = group_num;
                        else
                            group_num = group_num + 1;
                            dih_groups(d) = group_num;
                            start_dih = dihedrals(d,1:this_DOF);
                        end
                    end
                else
                    dih_groups = 1:size(dihedrals,1);
                end
                f = fopen(strcat(folder1, pdb_names, '_', Residue_names{to_write}, num2str(Residue_id(to_write)), '_dihedral_groups_', num2str(variant), '.txt'), 'w');
                for d = 1:size(dih_groups,1)
                    fprintf(f, '%f \n', dih_groups(d));
                end
                fclose(f);
                
            end
            f = fopen(strcat(folder1, pdb_names, '_parameters_', num2str(variant),'_g', num2str(group), '.txt'), 'w');
            if min_energy == 0
                fprintf(f, '%s %d %f \n', pdb_names, max_depth, 0.00001);
            else
                fprintf(f, '%s %d %f \n', pdb_names, max_depth, min_energy*1.01);
            end
            for res1 = 1:max_depth
                current_coord = eval(strcat('coord_', num2str(res1)));
                
                fprintf(f, '%s %d %f %d %d \n', Residue_names{res1}, Residue_id(res1), mins(res1),size(current_coord,1)-1, size(current_coord,2)-1 );
            end
            fclose(f);
            
            
            
            
            min_found = 1;
            
            else
                if min_tostart < 0.01
                 min_tostart = min_tostart*5;
                 min_found = 0;
                else
                    min_tostart = min_tostart*2;
                    min_found = 0;
                end
                
        end
        
        end
        end
end

toc
end
end


