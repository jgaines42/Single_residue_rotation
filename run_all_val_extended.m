%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = run_single_residue_in_dipeptide(PDB_name1, which_res,folder_name)
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

function [] = run_all_val_extended(folder_name, save_folder)

size_number = 1; %For SASA paper
next_pro = 0;
load('Val_coordinates.mat');
all_data = zeros(72,1);
for all_Val = 1:1000
    
    Position = Resi_coordinates(:,(all_Val-1)*3 + 1: all_Val*3);
    
    
    energy_cutoff = 1000000000;         %Initial energy cutoff
    
    %Set up the residue
    resiName = 'Val';
    resiId = 1;
    
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
    
    XtalPosition = Position;
    Atom_sizes = [1.88000000000000;1.61000000000000;1.42000000000000;1.64000000000000;1.88000000000000;1.61000000000000;1.42000000000000;1.88000000000000;1.88000000000000;1.88000000000000;0;0;0;0;0;0;0;0;0;1.64000000000000;1.88000000000000;0];
    
    
    
    
    %   sampledCombin = [sampledCombin, zeros(size(sampledCombin,1), DOF+2)];
    
    %% Get original chi and energy
    orig = [];
    
    res_name = resiName;
    switch_residue_setup;
    kbt = 1;
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
    %    save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_dipeptide_original_energy.mat'),'total_energy')
    %   save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId), '_dipeptide_original.mat'), 'orig')
    
    %%
    
    e_cutoff = .0001;
    energy = 10;
    
    %Run until a lowest energy value is lower than the cutoff
    
        [energy, dihedrals] = Find_all_energy_onlyChi1_dipeptide(Position, resiName, Atom_sizes,  0, e_cutoff);
        if sum(exp(-energy*10))>0
             all_data(:) = all_data(:)+exp(-energy*10)/sum(exp(-energy*10));
        end
    

end
            save(strcat(save_folder, PDB_name1, '_', resiName, num2str(resiId),'_dipeptide_single_rotation_minE.mat'), 'all_data');
            
        end
        



