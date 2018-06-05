%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [allDipeptide,next_pro, numAtoms, DOF] = isolate_dipeptide(tempModel2, all_ids, resiId resiName)
%
% Isolates a dipeptide in a protein. Returns the dipeptide, number of
% DOF and number of atoms in the dipeptide
%
% Input:
%   tempModel2: cell array of the entire PDB
%   all_ids: all residue ids of the protein
%   resId: the residue to be isolated
%   resiName: 3 letter abbreviation of the atom name
%
% Output:
%   allDipeptide: cell array of the dipeptide
%   next_pro: Indicates if the next residue is a Proline (which shortens
%   the length of the dipeptide). 1 = is proline
%
% Notes:
% A dipeptides is the amino acid being studied plus Calpha, C and O of the
% prior amino acid and N, H and Calpha of the next amino acid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [allDipeptide,next_pro, numAtom, DOF] = isolate_dipeptide(tempModel2, all_ids, resiId, resiName)

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


ind0=ismembc(all_ids,resiId);
tempAllIle = tempModel2(ind0, :);
ind0a = strcmp(tempAllIle(:,3), '');
ind0b = strcmp(tempAllIle(:,3), 'A');
ind0 = logical(ind0a+ind0b);
ind1 = strcmp(tempAllIle(:,3),'B');
model_B = tempAllIle(ind1,:);
tempAllIle = tempAllIle(ind0,:);


% itempFN = 4;
% find res -1, CA, C, O
ind0=ismembc(all_ids,resiId-1);
tempModel3=tempModel2(ind0,:);
ind0a = strcmp(tempModel3(:,2), 'CA');
ind0b = strcmp(tempModel3(:,2),'C');
ind0c  = strcmp(tempModel3(:,2),'O');

ind0 = logical(ind0a+ind0b+ind0c);
tempAllResm1=tempModel3(ind0,:);
ind0a = strcmp(tempAllResm1(:,3),'');
ind0b = strcmp(tempAllResm1(:,3),'A');
ind0 = logical(ind0a + ind0b);
tempAllResm1 = tempAllResm1(ind0,:);

% find res +1 N, CA, H
ind0=ismembc(all_ids,resiId+1);
tempModel3=tempModel2(ind0,:);
ind0a = strcmp(tempModel3(:,2), 'N');
ind0b = strcmp(tempModel3(:,2),'CA');
ind0c  = strcmp(tempModel3(:,2),'H');
ind0 = logical(ind0a+ind0b+ind0c);
tempAllResp1=tempModel3(ind0,:);

counter = 2;
while sum(ind0) ==0 && counter<max(resiId)
    ind0=ismembc(all_ids,resiId+counter);
    tempModel3=tempModel2(ind0,:);
    ind0a = strcmp(tempModel3(:,2), 'N');
    ind0b = strcmp(tempModel3(:,2),'CA');
    ind0c  = strcmp(tempModel3(:,2),'H');
    ind0 = logical(ind0a+ind0b+ind0c);
    tempAllResp1=tempModel3(ind0,:);
    counter = counter + 1;
end
next_pro = 0;
if sum(ind0) ~= 0
    %Deal with Pro
    if strcmp(tempAllResp1(1,4),'PRO')
        next_pro = 1;
    else
        next_pro = 0;
    end
    
    ind0a = strcmp(tempAllResp1(:,3),'');
    ind0b = strcmp(tempAllResp1(:,3),'A');
    ind0 = logical(ind0a + ind0b);
    tempAllResp1 = tempAllResp1(ind0,:);
end

allDipeptide=[tempAllResm1;tempAllIle;tempAllResp1];

end
