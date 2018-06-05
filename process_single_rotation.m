%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = process_single_rotation(what_res,resiName, save_folder, is_dipeptide)
% Processes the results of single rotation code (the version that returns
% the dihedrals with minimum energy and saves the distribution of dihedrals
% and the delta chi values. Details of the algorithm can be found in Gaines
% et al (2018) Proteins. 86, 581-591 and in the README file for this
% package
%
% Input:
% PDB_name: 4 residue PDB code
% what_res : Residue id
% resiName: 3 letter abbreviation of amino acid type
% save_folder: full path to folder with data, should end with '/'
% is_dipeptide: Whether or not code was run in context of dipeptide (1) or
% full protein (0)
%
% Output:
% *_all_dchi_values.mat file containing 50 delta chi values
% *_all_chi_distribution.mat file containing the frequency distribution for
% each of 50 combinations of run
%
% Notes:
% Currently only works on uncharged amino acids
% Requires the following files
% strcat(save_folder ,resiName, '_',  what_res, num2str( what_res), '_original.mat'));
% strcat(saved_folder, resiName, '_',what_res, num2str(what_res),'_dipeptide_single_rotation_minE.mat')
% or strcat(saved_folder, resiName, '_',what_res, num2str(what_res),'_protein_single_rotation_minE.mat')
% This file should contain the following columns:
% 1: which ba/bl varient
% 2: chi 1
% 3: chi 2
% continue to number of chi
% DOF+2 (last column): energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = process_single_rotation(PDB_name,what_res,resiName, save_folder, is_dipeptide)

KBT = 0.001; % Change this to make spheres squshier!
if is_dipeptide
    sub_name = '_dipeptide';
else
    sub_name = '_protein';
end
Resi_name = upper(resiName);

all_loc = [resiName,what_res];

num_babl_sampled = 300; % Total number of ba/bl that were run
num_samples = 50; % Number of delta chi values to produce
sample_size = 50; % Number of ba/bl runs to include in each analysis


rng(1); % Seed random number generator so results are constant


if exist(strcat(save_folder, PDB_name, '_', lower(resiName), num2str(what_res),sub_name,'_single_rotation_minE.mat'))
    load(strcat(save_folder, PDB_name, '_',lower(resiName), num2str(what_res),sub_name,'_single_rotation_minE.mat'));
    if size(unique(all_data(:,1)),1) == num_babl_sampled
        % If data existed and there are 300 unique runs, keep going
        
        % Convert Energy to probability  using Boltzman equation
        all_data(:,size(all_data,2)+1) = exp(-all_data(:,size(all_data,2))/KBT);
       
        % Determine how many degrees of freedom
        if ismember(Resi_name , {'ILE' ; 'LEU'; 'PHE'; 'TYR'; 'TRP'})
            DOF =2;
        elseif Resi_name == 'MET'
            DOF = 3;
        else
            DOF = 1;
        end
        
        best_dih = [];
        all_dchi = [];
        
        % Create num_samples groups of results
        for randloop = 1:num_samples
            to_sample = randperm(num_babl_sampled,sample_size);
            all_dih = [];
            
            % Extract data for each randomly chosen ba/bl variant add it to
            % existing data (with boltzman weighting)
            for j = 1:sample_size
                ind0 = ismember(all_data(:,1), to_sample(j));   % Find variant
                this_dih = all_data(ind0,2:size(all_data,2));   % Extract dihedral and energy
                if size(this_dih,1)>0
                    P_E = this_dih(1,size(this_dih,2));         % Extract boltzman weight
                    
                    % If this is the first variant, just add to all_dih
                    if j == 1
                        all_dih = this_dih(:,1:DOF);
                        all_dih(:,DOF+1) = P_E;
                    else % If not, see if dihedral is already in all_dih, if not add it.
                        
                        % This next if statement concatinates this_dih and all_dih and then finds unique rows.
                        % If the number of unique rows is the same as the initial number of items, than this_dih is not
                        % contained in all_dih and should be added to it.
                        if size(unique([this_dih(:,1:DOF);all_dih(:,1:DOF)], 'rows'),1) == (size(this_dih,1) + size(all_dih,1))
                            all_dih = [all_dih;this_dih(:,1:DOF), repmat(P_E,size(this_dih,1),1)];
                            
                        else
                            % Determine which members of this_dih are in all_dih and add their P_E to all_dih
                            ind0 = ismember(all_dih(:,1:DOF), this_dih(:,1:DOF),'rows');
                            all_dih(ind0,DOF+1) = all_dih(ind0,DOF+1)+P_E;
                            % Concatinate the rest of this_dih to all_dih
                            ind1 = ismember(this_dih(:,1:DOF), all_dih(:,1:DOF), 'rows');
                            all_dih = [all_dih; this_dih(~ind1,1:DOF), repmat(P_E, sum(~ind1),1)];
                        end
                    end
                end
            end
            
            % Once we have sampled all sample_size ba/bl variants, we find
            % the most probable (weighted) dihedral and calculated delta chi
            max_val = max(all_dih(:,DOF+1));
            ind0 = ismember(all_dih(:,DOF+1), max_val); % Find all dihedral with that value
            best_dih = [best_dih;all_dih(ind0,1:DOF),repmat(randloop,sum(ind0),1),  find(ind0 == 1)];
            
            % If there are multiple dihedral angles with the same
            % probability, randomly choose one to use
            chi1 = all_dih(ind0,1:DOF);
            rand_numb1 = randi(size(chi1,1), 1);
            chi1 = chi1(rand_numb1,:);
            
            %Calculate delta chi
            this_res = calc_delta_chi(PDB_name, resiName, what_res, save_folder,DOF, chi1);

            all_dchi = [all_dchi;this_res];
        end
        save(strcat(save_folder,PDB_name, '_', lower(resiName),num2str(what_res), sub_name,'_best_dih.mat'), 'best_dih');

        save(strcat(save_folder,PDB_name, '_', lower(resiName),num2str(what_res), sub_name,'_all_dchi_values.mat'), 'all_dchi');
    else
        error('300 variants were not run')
    end
   
else
    error('Input files do not exist')
end


end