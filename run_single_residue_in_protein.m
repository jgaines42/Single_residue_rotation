%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = run_single_residue_in_protein(PDB_name1, which_res,folder_name)
% Performs rotations on a single residue. Save the lowest energy state of
% 100 differen bond length and angle variants
%
% Input:
%   PDB_name1: the name of the PDB to be run. Should be saved as
%       XXX.mat where PDB_name1 = XXX
%   which_res: the residue ID
%   folder_name: The full path name to the folder containing the PDB file.
%   save_folder: Folder to save results to
% Output:
%   XXX_original.mat: the original chi values of the residue
%   XXX_original_energy.mat: the original energy of the residue
%   XXX_single_rotation_minE.mat: lowest energy state(s) for each ba/bl
%       variant
%       column 1: variant
%       column 2: minimum energy
%       column 3-n: chi values at that energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = run_single_residue_in_protein(PDB_name1, which_res,folder_name, save_folder)

size_number = 71; %For SASA paper

% Load PDB
PDB_name = strcat(folder_name, lower(PDB_name1), '.mat');
load(PDB_name);
ind1 = strcmp(tempModel2(:,3),'B');
tempModel2 = tempModel2(~ind1,:);
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);

res_ids = cell2mat(tempModel2(:,6));

%{
%relabel residue ids of second chain so all > ids of 1st chain
for i = 2:size(tempModel2,1)
    if res_ids(i,1) < res_ids(i-1,1)
        res_ids(i:size(tempModel2,1),1) = res_ids(i:size(tempModel2,1),1) + 1000;
    end
end
tempModel2(:,6) = num2cell(res_ids);
%}
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);

%Add sizes to everything
tempModel2 = add_sizes_protein(tempModel2,size_number);
sizes_all = cell2mat(tempModel2(:,12));


energy_cutoff = 1000000000;         %Initial energy cutoff

%Set up the residue
ind0 = find(res_ids == which_res);
resiName = tempModel2{ind0(1),4}; %Get name of interface resiude
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


%Set up the dihedral angle
%sampledCombin = set_up_dihedral_angles(DOF);

%Isolate dipeptide
[allDipeptide,next_pro] = isolate_dipeptide(tempModel2, res_ids, resiId);


if DOF >0 && (size(allDipeptide,1) == numAtom || (size(allDipeptide,1) == numAtom-1 && next_pro ==1))
    
    %Check the order of the dipeptide
    if next_pro == 0
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-3,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-3,:)= new_Dipeptide;
    else
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-2,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-2,:)= new_Dipeptide;
    end
    ind0 = ismember(cell2mat(tempModel2(:,1)), cell2mat(allDipeptide(:,1)));
    rest_of_pro = tempModel2(~ind0,:);
    ind0 = ismember(cell2mat(rest_of_pro(:,6)), resiId);
    rest_of_pro = rest_of_pro(~ind0,:);
    ind0 = ismember(rest_of_pro(:,3), {'','A'});
    rest_of_pro = rest_of_pro(ind0,:);
    
    
    %% Save dipeptide positions
    XtalPosition = cell2mat(allDipeptide(:,8:10));
    Position = XtalPosition;
    Atom_sizes = cell2mat(allDipeptide(:,12));
    
    
    %% Limit Rest of Pro to only be atoms within 10A of Cb\
    % basic checks have shown that Phe, Met and Leu always have all
    % atoms stay within 6A of Cb, so add 4 just in case
    pro_pos = cell2mat(rest_of_pro(:,8:10));
    distp = repmat(Position(8,:),size(pro_pos,1),1)-pro_pos;
    distempP = sqrt(sum(distp.^2,2));
    ind1 = distempP <= 20;
    rest_of_pro = rest_of_pro(ind1,:);
    
    if correct_now == 1 && DOF >0
     %   sampledCombin = [sampledCombin, zeros(size(sampledCombin,1), DOF+2)];
        
        %% Get original chi and energy
        orig = [];
        
        rest_pro_position = cell2mat(rest_of_pro(:,8:10));
        res_name = resiName;
        switch_residue_setup;
        kbt = 1;
        set_up_clashArrays;
        total_energy = 0;
        if DOF >= 1
            orig = [InitChi1];
            [chi1_energy] = get_energy_wProtein(c1_Clash, Position,  protein_clash_c1, rest_pro_position);
            total_energy = total_energy +chi1_energy;
            
            if DOF >= 2
                orig = [InitChi1, InitChi2];
                [chi2_energy] = get_energy_wProtein(c2_Clash, Position,  protein_clash_c2, rest_pro_position);
                total_energy = total_energy + chi2_energy;
            end
            if DOF >= 3
                orig = [InitChi1, InitChi2, InitChi3];
                [chi3_energy] = get_energy_wProtein(c3_Clash, Position,  protein_clash_c3, rest_pro_position);
                total_energy = total_energy + chi3_energy;
            end
        end
        if CH3 == 2
            [HG2_energy] = get_energy_wProtein(HG2_Clash, Position,  protein_clash_HG2, rest_pro_position);
            total_energy = total_energy + HG2_energy;
        end
        if CH3 >=1
            [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  protein_clash_HG1, rest_pro_position);
            total_energy = total_energy + HG1_energy;
        end
        if OH ==1
            [OH_energy] = get_energy_wProtein(OH_Clash, Position,  protein_clash_OH, rest_pro_position);
            total_energy = total_energy + OH_energy;
        end
        save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_original_energy.mat'),'total_energy')
        save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_original.mat'), 'orig')
        
        %%
        %Get bond length and angle variants
        [Position_100, rest_of_pro, total_found] = get_1000(allDipeptide, resiName,rest_of_pro, next_pro);
        if total_found >= 300
            all_data = [];
            
            %Run 100 variants
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
                    [energy, dihedrals] = Find_lowest_energy_1Residue_varyBABL(Position, resiName, Atom_sizes,  rest_of_pro, next_pro, e_cutoff);
                end
                if energy < e_cutoff
                    all_data = [all_data; repmat([variations, energy],size(dihedrals,1),1), dihedrals];
                end
                
            end
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId),'_single_rotation_minE.mat'), 'all_data');
            
        end
        
    end
    
end

end


