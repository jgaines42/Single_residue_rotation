%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [energy] = get_energy_wProtein(DA_Clash, Position,  protein_clash, Pro_Position)
%
% Calculates the energy between a set of atoms and the rest of the protein
% as indicated in DA_Clash and protein_clash
%
% Input:
%   DA_Clash: Indexes of the atom pairs between the moving side chain and
%       the rest of the Dipeptidde. Indexed based on the Dipeptide
%   Position: Position of the Dipeptide
%   protein_clash: Indexes of the atom pairs between the moving side chain
%       and the rest of the protein. Column 1 is indexed based on the
%       Dipeptide and Column 2 is indexed based on the protein
%   Pro_Position: Position of the rest of the protein
%
%  Output:
%       energy: RLJ energy due to overlaps
%
%  Note: 
%   Equation use is Repulsive Lennard-Jones euqation. The inner power in the
%   equation is 3 because DA_Clash(:,3) and distemp are sigma_ij^2 and
%   r_ij^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [energy] = get_energy_wProtein(DA_Clash, Position,  protein_clash, Pro_Position)

% Calculate index of clashes within the dipeptide
if size(DA_Clash,1) >0
    dist = Position(DA_Clash(:,1),:)-Position(DA_Clash(:,2),:);
    distemp = sum(dist.^2,2);
    clash = distemp - DA_Clash(:,3);
    ind_clash = find(clash < 0);
else
    ind_clash = [];
end

%Calculate index of clashes of the dipeptide with the protein
if size(protein_clash,1)>0
    distp = Position(protein_clash(:,1),:)-Pro_Position(protein_clash(:,2),:);
    distempP = sum(distp.^2,2);
    clashP = distempP - protein_clash(:,3);
    ind_clashP = find(clashP < 0);
else
    ind_clashP = [];
end

%Calculate energy
if size(DA_Clash,1) >0 && size(protein_clash,1)>0
    if size(ind_clash,1) > 0 || size(ind_clashP,1) >0
        energy = sum((1 - (DA_Clash(ind_clash,3)./distemp(ind_clash)).^3).^2) + sum((1 - (protein_clash(ind_clashP,3)./distempP(ind_clashP)).^3).^2);
    else
        energy = 0;
    end
elseif size(protein_clash,1)>0 % If just doing clash with protein, not peptide
    if  size(ind_clashP,1) >0
        energy =   sum((1 - (protein_clash(ind_clashP,3)./distempP(ind_clashP)).^3).^2);
    else
        energy = 0;
    end
elseif size(DA_Clash,1) >0
    if size(ind_clash,1) > 0
        energy = sum((1 - (DA_Clash(ind_clash,3)./distemp(ind_clash)).^3).^2) ;
    else
        energy = 0;
    end
end
end