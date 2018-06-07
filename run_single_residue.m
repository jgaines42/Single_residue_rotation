%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = run_single_residue(PDB_name1, which_res,folder_name, save_folder, is_dipeptide)
% Performs rotations on a single residue. Save the lowest energy state of
% 300 differen bond length and angle variants
%
% Input:
%   PDB_name1: the name of the PDB to be run. Should be saved as
%       XXXX_H.pdb or XXXX.mat where PDB_name1 = XXXX
%   which_res: the residue ID
%   folder_name: The full path name to the folder containing the PDB file.
%   save_folder: Folder to save results to
%   is_dipeptide: 1 to run just the dipeptide, 0 to run in context of the
%   protein
% Output:
%   X_dipeptide_original.mat: the original chi values of the residue (or
%       _protien_)
%   X_dipeptide_original_energy.mat: the original energy of the residue (or
%       _protien_)
%   X_dipeptide_single_rotation_minE.mat: lowest energy state(s) for each ba/bl
%       variant (or _protien_)
%       column 1: variant
%       column 2 to (n-1): chi values at minimum energy
%       column n: minimum energy
%   X_dipeptide_all_dchi_values.mat: The dchi values for samplings of the lowest
%   energy values found (see README for algorithm details) (or _protien_)
% Notes:
% Protein must be saved in file XXXX_H.pdb created
% using download_preprocess_pdb.py which also adds the hydrogen atoms to the protein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = run_single_residue(PDB_name1, which_res,folder_name, save_folder,is_dipeptide)

size_number = 7; %Used in add_sizes_protein().

% Load PDB.mat file if it exists, if not, create it
if ~exist(strcat(folder_name, PDB_name1, '.mat'))
    if ~exist(strcat(folder_name, PDB_name1, '_H.pdb'))
        error('PDB file with hydogen atoms added must exist')
    else
        pdbstruct = pdbread(strcat(folder_name, PDB_name1, '_H.pdb'));
        x = pdbstruct.Model.Atom;
        tempModel1=struct2cell(x);
        tempModel2=reshape(tempModel1,size(tempModel1,1),size(tempModel1,3))';
        save(strcat(folder_name, PDB_name1,'.mat'), 'tempModel2');
    end
else
    load(strcat(folder_name, PDB_name1,'.mat'));
end

%Renumber atoms to be sequential
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);

%relabel residue ids of second chain so all > ids of 1st chain
res_ids = cell2mat(tempModel2(:,6));
for i = 2:size(tempModel2,1)
    if res_ids(i,1) < res_ids(i-1,1)
        res_ids(i:size(tempModel2,1),1) = res_ids(i:size(tempModel2,1),1) + 1000;
    end
end
tempModel2(:,6) = num2cell(res_ids);
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);

%Add sizes to everything
tempModel2 = add_sizes_protein(tempModel2,size_number);
sizes_all = cell2mat(tempModel2(:,12));

%Set up the residue
ind0 = find(res_ids == which_res);
resiName = tempModel2{ind0(1),4}; %Get name of amino acid
resiName(2:3) = lower(resiName(2:3));
resiId = double(which_res);


%Isolate dipeptide
[allDipeptide,next_pro, numAtom, DOF] = isolate_dipeptide(tempModel2, res_ids, resiId, resiName);

% Make sure the full resiue was present
if DOF >0 && (size(allDipeptide,1) == numAtom || (size(allDipeptide,1) == numAtom-1 && next_pro ==1))
    
    %Put the atoms in the correct order
    if next_pro == 0
        [new_Dipeptide, correct_now, ~] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-3,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-3,:)= new_Dipeptide;
    else
        [new_Dipeptide, correct_now, ~] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-2,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-2,:)= new_Dipeptide;
    end    
    
    %% Save dipeptide positions
    XtalPosition = cell2mat(allDipeptide(:,8:10));
    Position = XtalPosition;
    Atom_sizes = cell2mat(allDipeptide(:,12));
    
    %% get rest of the protein
    if is_dipeptide
        rest_of_pro = [];
        rest_of_pro_position = [];
    else
        ind0 = ismember(cell2mat(tempModel2(:,1)), cell2mat(allDipeptide(:,1)));
        rest_of_pro = tempModel2(~ind0,:);
        ind0 = ismember(cell2mat(rest_of_pro(:,6)), resiId);
        rest_of_pro = rest_of_pro(~ind0,:);
        ind0 = ismember(rest_of_pro(:,3), {'','A'});
        rest_of_pro = rest_of_pro(ind0,:);
        
        %% Limit Rest of Pro to only be atoms within 10A of Cb\
        % basic checks have shown that Phe, Met and Leu always have all
        % atoms stay within 6A of Cb, so add 4 just in case
        pro_pos = cell2mat(rest_of_pro(:,8:10));
        distp = repmat(Position(8,:),size(pro_pos,1),1)-pro_pos;
        distempP = sqrt(sum(distp.^2,2));
        ind1 = distempP <= 20;
        rest_of_pro = rest_of_pro(ind1,:);
        rest_of_pro_position = cell2mat(rest_of_pro(:,8:10));

    end
   
    %If the correct atom ordering was possible
    if correct_now == 1
        
        %% Finish setup and get original chi and energy
         Amino_acid = switch_residue_setup(resiName);
         if is_dipeptide
             
            [Init_dihedrals, Clash_lists] = set_up_clashArrays_dipeptide( resiName,next_pro,Atom_sizes, Amino_acid, Position);
         else
             [Init_dihedrals, Clash_lists] = set_up_clashArrays_protein( resiName,next_pro,Atom_sizes, Amino_acid, Position,rest_of_pro);
         end
         [total_energy] = get_original_energy(Clash_lists, Position, Amino_acid.DOF, Amino_acid.CH3, Amino_acid.OH,rest_of_pro_position);
       
        if DOF == 1
            orig = [Init_dihedrals.InitChi1];
        elseif DOF == 2
            orig = [Init_dihedrals.InitChi1,Init_dihedrals.InitChi2];
        elseif DOF == 3
            orig = [Init_dihedrals.InitChi1,Init_dihedrals.InitChi2,Init_dihedrals.InitChi3];
        end
        if is_dipeptide
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_dipeptide_original_energy.mat'),'total_energy')
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_dipeptide_original.mat'), 'orig')
        else
            
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_protein_original_energy.mat'),'total_energy')
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_protein_original.mat'), 'orig')
        end
        %% Get bond length and angle variants

        if is_dipeptide %create empty array for get_1000 to work
            temp_rest = cell(1,15); % this line and the next allows get_1000 to work for the dipeptide model
            temp_rest(1,8:10) = num2cell([0 0 0]);
        else
            temp_rest = rest_of_pro;
        end
        
        [Position_100, rest_of_pro, total_found] = get_1000(allDipeptide, resiName,temp_rest, next_pro, Amino_acid);
        if is_dipeptide
            rest_of_pro_position = [];
        else
            rest_of_pro_position = cell2mat(rest_of_pro(:,8:10));
        end
        %% If we were able to get the 300 variants we need, keep runing
        if total_found >= 300
            all_data = [];

            %Run 300 variants
            for variations = 1:300
                Position = Position_100(:,(variations-1)*3 + 1: variations*3); %Extract variant
                energy = 100;   %Set initial energy to be high
                
                % Calculate new initial dihedral angles
                if DOF == 3
                    InitChi3=calcDA2( Amino_acid.iChi3Array,Position); %%%%%%%%
                    Init_dihedrals.InitChi3=mod(real(InitChi3),360);
                end
                if DOF >= 2
                    
                    InitChi2=calcDA2( Amino_acid.iChi2Array,Position); %%%%%%%%
                    Init_dihedrals.InitChi2=mod(real(InitChi2),360);
                end
                if DOF >= 1
                    InitChi1=calcDA2( Amino_acid.iChi1Array,Position);
                    Init_dihedrals.InitChi1=mod(real(InitChi1),360);
                end
                if Amino_acid.CH3 == 2
                    InitChi_HG2=calcDA2( Amino_acid.HG_Array_2,Position);
                    Init_dihedrals.InitChi_HG2=mod(real(InitChi_HG2),360);
                end
                if Amino_acid.CH3 >= 1
                    InitChi_HG1=calcDA2( Amino_acid.HG_Array_1,Position);
                    Init_dihedrals.InitChi_HG1=mod(real(InitChi_HG1),360);
                end
                if Amino_acid.OH == 1
                    InitOH=calcDA2( Amino_acid.iOHArray,Position); %InitChi1
                    Init_dihedrals.InitOH=mod(real(InitOH),360);
                end

                %Set initial cutoff based off of old runs
                if variations >1
                    if max(all_data(:,DOF+2) >0)
                        e_cutoff = max(all_data(:,DOF+2))/100;
                    else
                        e_cutoff = 10^-5;
                    end
                else
                    e_cutoff = .0001;
                end
                %Run until a lowest energy value is lower than the cutoff
                while energy >= e_cutoff
                    e_cutoff = e_cutoff *10;
                    
                    [energy, dihedrals] = Find_lowest_energy_1Residue_varyBABL(Position, e_cutoff, Clash_lists, Amino_acid, Init_dihedrals, rest_of_pro_position);
                    
                end
                if energy < e_cutoff
                    all_data = [all_data; repmat([variations],size(dihedrals,1),1) dihedrals,repmat([ energy],size(dihedrals,1),1)];
                end
                
            end
            if is_dipeptide
                save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId),'_dipeptide_single_rotation_minE.mat'), 'all_data');
            else
                save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId),'_protein_single_rotation_minE.mat'), 'all_data');
            end 
            %Process the data
            process_single_rotation(PDB_name1,resiId,resiName, save_folder, is_dipeptide);
        else
            error('Did not find enough variants')
        end
    else
        error('Atoms are not labeled correctly')
    end
else
    error('Residue is missing atoms')
end

end


