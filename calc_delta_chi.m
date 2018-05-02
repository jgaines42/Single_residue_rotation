%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function this_res = calc_delta_chi(PDB_name, Resi_name, Resi_id, save_folder, DOF, chi_values)
% Calcualtes delta chi for the given amino acid 
%
% Input:
% PDB_name: 4 letter PDB code
% Resi_name: 3 letter amino acid name (in caps)
% Resi_id: residue id
% save_folder: full path to the folder witht the data with / at end
% DOF: number of side chain dihedral angles
% chi_values: The dihedral angle(s) value to compare to the crystal structure
%
% Output:
% this_res: Delta chi for chi_values
%
% Notes:
% Requires crystal structure chi values to be stored in *_original.mat
% Only considers chi2 from 0 to 180 for Phe and Tyr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function this_res = calc_delta_chi(PDB_name, Resi_name, Resi_id, save_folder, DOF, chi_values)
    load(strcat(save_folder ,PDB_name, '_',  Resi_name, num2str( Resi_id), '_original.mat'));
   
    if DOF == 2
        % If Phe or Tyr, we keep dihedral within 0 and 180
        if (Resi_name == 'PHE' | Resi_name== 'TYR') & orig(2) > 180
            orig(2) = orig(2)-180;
        end
        d1 = abs(chi_values(:,1) - orig(1));
        d2 = abs(chi_values(:,2) - orig(2));
        ind0 = d1 > 180;
        d1(ind0) = 360-d1(ind0);
        ind0 = d2 > 180;
        d2(ind0)  = 360-d2(ind0);
        this_res = min(sqrt(d1.^2 + d2.^2));

    elseif DOF == 1

        d1 = abs(chi_values(1)-orig(1));
        ind0 = d1 > 180;
        d1(ind0) = 360-d1(ind0);
        this_res = min(d1);


    elseif DOF == 3
        d1 = abs(chi_values(:,1) - orig(1));
        d2 = abs(chi_values(:,2) - orig(2));
        d3 = abs(chi_values(:,3) -orig(3));
        ind0 = d1 > 180;
        d1(ind0) = 360-d1(ind0);
        ind0 = d2 > 180;
        d2(ind0)  = 360-d2(ind0);
        ind0 = d3 > 180;
        d3(ind0)  = 360-d3(ind0);
        this_res = min(sqrt(d1.^2 + d2.^2 + d3.^2));


    end
end