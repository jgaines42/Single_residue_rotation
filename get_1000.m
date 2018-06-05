%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Position_1000, rest_of_pro, total_found] =get_1000(allDipeptide, resiName, rest_of_pro, next_pro, Amino_acid)
%   
% Generates up to 1000 instances of a given residue, taken from Dunbrack
% 1.7.
%
% Input: 
%   allDipeptide: Cell array of the dipeptide
%   resiName: Residue name
%   rest_of_protein: The rest of the protein (minus the dipeptide)
%   next_pro: 1 if next residue is Proline
%   Amino_acid: Structure of data of the amino acid type
%
% Output:
%   Position_1000: coordiant array of 1000 instances of the residue. Each
%       set of 3 columns is an x,y,z pair
%   rest_of_pro: New rest_of_pro cell array (because everything gets moved)
%   total_found: The number of instances of each residue that were
%       extracted
%
% Notes:
%   The coordinate array loaded has removed all instances of the residue
%   that have bond lengths or angles outside of 3 standard deviations of
%   the original distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Position_1000, rest_of_pro, total_found] =get_1000(allDipeptide, resiName, rest_of_pro, next_pro, Amino_acid)

Position_1000 = zeros(size(allDipeptide,1), 1000*3);
residues = 1;

bond_angles_all = [];
bond_lengths_all = [];
location = [1]; % Stores where to start looking for the next coordinates
Model_new = cell2mat(rest_of_pro(:,8:10));

%%
load(strcat(resiName, '_coordinates.mat')) % This file contains the coordinates of other instances of the amino acid

numAtom = Amino_acid.numAtom;
iChi1Array = Amino_acid.iChi1Array;
moveAtomID2 = Amino_acid.moveAtomID2;

if next_pro == 1
    numAtom = numAtom-1;
end
%% Save dipeptide coordinates
XtalPosition = cell2mat(allDipeptide(:,8:10));
Position = XtalPosition;
Position_1000(:,1:3) = Position;

%Rotate rest of protein in same way that dipeptide will be rotated
TempPosition=Position-repmat(Position(5,:),numAtom,1);
Model_new = Model_new - repmat(Position(5,:), size(Model_new,1), 1);

% Calculate what dihedral C-N-Ca-Cb angle should be
InitAng=calcDA2([2 4 5 8],TempPosition);
InitAng=mod(real(InitAng),360);

%Calculate Chi1
InitChi1=calcDA2(iChi1Array,TempPosition);
InitChi1=mod(real(InitChi1),360);


%% Orient original coordinates
%Rotate everything around y axis until N is on x-y plane
rot_Ny =  atan2d(TempPosition(4,3), TempPosition(4,1));
R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
TempPosition = (R*TempPosition')';
Model_new = (R*Model_new')';


%Rotate everything around z axis until N is on x axis
rot_Ny =  -atan2d(TempPosition(4,2), TempPosition(4,1));
R =    [ cosd(rot_Ny) ,-sind(rot_Ny),0;sind(rot_Ny), cosd(rot_Ny),0;0 0 1];
TempPosition = (R*TempPosition')';
Model_new = (R*Model_new')';
if TempPosition(4,1) < 0
    rot_Ny =  180;
    R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
    TempPosition = (R*TempPosition')';
    Model_new = (R*Model_new')';
end

%Now N is on the x axis, move Cb to the x,y plane, rotate on x
rot_Ny =  -atan2d(TempPosition(8,3), TempPosition(8,2));
R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
TempPosition = (R*TempPosition')';
Model_new = (R*Model_new')';
if TempPosition(8,2) < 0
    rot_Ny = 180;
    R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
    TempPosition = (R*TempPosition')';
    Model_new = (R*Model_new')';
end


rest_of_pro(:,8:10) = num2cell(Model_new);
Position_1000(:,1:3) = TempPosition;

TempPostion_orig = TempPosition;

%% Loop and get all new coordinates
i = 1;
for multi_pdb = 2:1000 %Make dipeptides with new bl/ba
    
    TempPosition = TempPostion_orig;
    
    %Grab new coordinates
    new_res = Resi_coordinates(:,[(i-1)*3+1 : i*3]);
    if next_pro ==1
        new_res = new_res(1:size(new_res,1)-1, :);
    end
 
   %% Orient new coordinates
    % Move Ca to zero
    new_res = new_res-repmat(new_res(5,:), numAtom,1);
    %Rotate everything around y axis until N is on x-y plane
    rot_Ny =  atan2d(new_res(4,3), new_res(4,1));
    R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
    new_res = (R*new_res')';
    
    %Rotate everything around z axis until N is on x axis
    rot_Ny =  -atan2d(new_res(4,2), new_res(4,1));
    R =    [ cosd(rot_Ny) ,-sind(rot_Ny),0;sind(rot_Ny), cosd(rot_Ny),0;0 0 1];
    new_res = (R*new_res')';
    
    %Check that N is on positive x axis
    if new_res(4,1) < 0
        rot_Ny =  180;
        R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
        new_res = (R*new_res')';
    end
    
    %Rotate everything so that Cb is in positive y, xy plane
    rot_Ny =  -atan2d(new_res(8,3), new_res(8,2));
    R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
    new_res = (R*new_res')';
    if new_res(8,2) < 0
        rot_Ny = 180;
        R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
        new_res = (R*new_res')';
    end
     
    
    %Move new sidechain coordinates to original backbone
    TempPosition([8,moveAtomID2],:) = new_res([8, moveAtomID2],:);
    
    %Rotate to the correct dihedral angle
    setAng = InitAng;
    newVal_ang = calcDA2([2 4 5 8],TempPosition);
    newVal_ang = mod(real(newVal_ang),360);
    
    subtract_array_2 = repmat(TempPosition(4,:),numAtom,1);
    delta_term_1 = pi*sign(newVal_ang)*newVal_ang/180;
    TempPosition = Rotate_DA(TempPosition, setAng, subtract_array_2, delta_term_1, [2 4 5 8], [8, moveAtomID2]);
    
    
    current_Chi1 = calcDA2(iChi1Array, TempPosition);
    current_Chi1 = mod(real(current_Chi1), 360);
    subtract_array_1 = repmat(TempPosition(iChi1Array(2),:),numAtom,1);
    delta_term_1 =  pi*sign(current_Chi1)*current_Chi1/180;
    TempPosition = Rotate_DA(TempPosition, InitChi1, subtract_array_1, delta_term_1, iChi1Array, moveAtomID2);
    
    Position_1000(:,(multi_pdb-1)*3+1:multi_pdb*3) = TempPosition;
    
    i = i+1;
end
total_found = multi_pdb;
end