%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [allDipeptide,next_pro] = isolate_dipeptide(tempModel2, all_ids, resiId)
%
% Isolates a dipeptide in a protein. Returns the dipeptide 
%
% Input:
%   tempModel2: cell array of the entire PDB
%   all_ids: all residue ids of the protein
%   resId: the residue to be isolated
%
% Output:
%   allDipeptide: cell array of the dipeptide
%   next_pro: Indicates if the next residue is a Proline (which shortens
%   the length of the dipeptide). 1 = is proline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [allDipeptide,next_pro] = isolate_dipeptide(tempModel2, all_ids, resiId)

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
