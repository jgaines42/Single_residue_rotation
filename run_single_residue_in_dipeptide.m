%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = run_single_residue_in_dipeptide(PDB_name1, which_res,folder_name)
% Performs rotations on a single residue. Save the lowest energy state of
% 100 differen bond length and angle variants
%
% Input:
%   PDB_name1: the name of the PDB to be run. Should be saved as
%       XXXX.mat where PDB_name1 = XXXX
%   which_res: the residue ID
%   folder_name: The full path name to the folder containing the PDB file.
%   save_folder: Folder to save results to
% Output:
%   X_original.mat: the original chi values of the residue
%   X_original_energy.mat: the original energy of the residue
%   X_single_rotation_minE.mat: lowest energy state(s) for each ba/bl
%       variant
%       column 1: variant
%       column 2: minimum energy
%       column 3-n: chi values at that energy
%   X_all_dchi_values.mat: The dchi values for samplings of the lowest
%   energy values found (see README for algorithm details
% Notes:
% Protein must be saved in file XXXX_H.pdb created
% using download_preprocess_pdb.py which also adds the hydrogen atoms to the protein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = run_single_residue_in_dipeptide(PDB_name1, which_res,folder_name, save_folder)

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

switch (resiName)
    case 'Ala'
        numAtom = 16;
        DOF = 0;
    case 'Ile'
        numAtom =25;
        DOF = 2;
        
    case 'Leu'
        numAtom =25;
        DOF = 2;
        
    case 'Val'
        numAtom = 22;
        DOF = 1;
        
    case 'Phe'
        numAtom = 26;
        DOF =2;
        
    case 'Trp'
        numAtom = 30;
        DOF = 2;
    case 'Tyr'
        numAtom = 27;
        DOF = 2;
    case 'Asn'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'Cys'
        DOF = 1;
        numAtom = 17;
        
    case 'Glu'
        DOF = 3;
        numAtom = 21;
        DOF = 0;
    case 'Met'
        numAtom = 23;
        DOF = 3;
        
    case 'Ser'
        numAtom = 17;
        DOF = 1;
    case 'Thr'
        numAtom = 20;
        DOF = 1;
    case 'Asp'
        numAtom = 18;
        DOF = 2;
        DOF = 0;
        
    case 'Gln'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'Arg'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'His'
        numAtom = 22;
        DOF = 2;
    case 'Lys'
        numAtom = 28;
        DOF = 4;
        DOF = 0;
    case 'Gly'
        DOF = 0;
        fprintf('Not yet supported\n' );
    case 'Pro'
        fprintf('Not yet supported\n' );
        DOF = 0;
    otherwise
        fprintf('Invalid amino acid\n' );
        DOF = 0;
end

%Isolate dipeptide
[allDipeptide,next_pro] = isolate_dipeptide(tempModel2, res_ids, resiId);

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
    
   
    %If the correct atom ordering was possible
    if correct_now == 1
        
        %% Get original chi and energy
        orig = [];
        
        res_name = resiName;
        switch_residue_setup;
        set_up_clashArrays_dipeptide;
        total_energy = 0;
        if DOF >= 1
            orig = [InitChi1];
            [chi1_energy] = get_energy_wProtein(c1_Clash, Position,  [], []);
            total_energy = total_energy +chi1_energy;
            
            if DOF >= 2
                orig = [InitChi1, InitChi2];
                [chi2_energy] = get_energy_wProtein(c2_Clash, Position,  [], []);
                total_energy = total_energy + chi2_energy;
            end
            if DOF >= 3
                orig = [InitChi1, InitChi2, InitChi3];
                [chi3_energy] = get_energy_wProtein(c3_Clash, Position,  [], []);
                total_energy = total_energy + chi3_energy;
            end
        end
        if CH3 == 2
            [HG2_energy] = get_energy_wProtein(HG2_Clash, Position,  [], []);
            total_energy = total_energy + HG2_energy;
        end
        if CH3 >=1
            [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  [], []);
            total_energy = total_energy + HG1_energy;
        end
        if OH ==1
            [OH_energy] = get_energy_wProtein(OH_Clash, Position,  [],[]);
            total_energy = total_energy + OH_energy;
        end
        save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_original_energy.mat'),'total_energy')
        save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_original.mat'), 'orig')
        
        %%
        temp_rest = cell(1,15);
        temp_rest(1,8:10) = num2cell([0 0 0]);
        %Get bond length and angle variants
        [Position_100, rest_of_pro, total_found] = get_1000(allDipeptide, resiName,temp_rest, next_pro);
        if total_found >= 300
            all_data = [];
            
            %Run 300 variants
            for variations = 1:300
                Position = Position_100(:,(variations-1)*3 + 1: variations*3);
                
                energy = 100;
                
                %Set initial cutoff based off of old runs
                if variations >1
                    if max(all_data(:,2) >0)
                        e_cutoff = max(all_data(:,2))/100;
                    else
                        e_cutoff = 10^-5;
                    end
                else
                    e_cutoff = .0001;
                end
                
                %Run until a lowest energy value is lower than the cutoff
                while energy >= e_cutoff
                    e_cutoff = e_cutoff *10;
                    [energy, dihedrals] = Find_lowest_energy_1Residue_varyBABL_dipeptide(Position, resiName, Atom_sizes,   next_pro, e_cutoff);
                end
                if energy < e_cutoff
                    all_data = [all_data; repmat([variations],size(dihedrals,1),1) dihedrals,repmat([ energy],size(dihedrals,1),1)];
                end
                
            end
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId),'_dipeptide_single_rotation_minE.mat'), 'all_data');
            
            %Process the data
            process_single_rotation(PDB_name1,resiId,resiName, save_folder, 1);
        end
    else
        error('Atoms are not labeled correctly')
    end
else
    error('Residue is missing atoms')
end

end


