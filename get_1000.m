%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Position_1000, rest_of_pro, total_found] =get_1000(allDipeptide, resiName, rest_of_pro, next_pro)
%   
% Generates up to 1000 instances of a given residue, taken from Dunbrack
% 1.7.
%
% Input: 
%   allDipeptide: Cell array of the dipeptide
%   resiName: Residue name
%   rest_of_protein: The rest of the protein (minus the dipeptide)
%   next_pro: 1 if next residue is Proline
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



function [Position_1000, rest_of_pro, total_found] =get_1000(allDipeptide, resiName, rest_of_pro, next_pro)

Position_1000 = zeros(size(allDipeptide,1), 1000*3);
residues = 1;

%pl_mi = 15; % This is how different the Phi/Phi of the new coordinates can be from the old ones

bond_angles_all = [];
bond_lengths_all = [];
location = [1]; % Stores where to start looking for the next coordinates
Model_new = cell2mat(rest_of_pro(:,8:10));

%%
load(strcat(resiName, '_coordinates.mat')) % This file contains the coordinates of other instances of the amino acid

switch (resiName)
    % The switch statement creates the following variables
    % Phi_ind = index of atoms that create Phi within the dipeptide
    % Psi_ind = index of atoms that create Psi within the dipeptide
    % iChi1Array = index of atoms that create Chi1 in dipeptide
    % moveAtomID2 = index of atoms that move when Chi1 changes
    
    case 'Ala'
        numAtom = 16;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 14];
        iChi1Array = [4 5 8 11];
        moveAtomID2 = [11:13];
    case 'Ile'
        numAtom =25;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 23];
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,11];
        moveAtomID=[11,15:16,20:22]; %%%%%%%% rotate chi2
        moveAtomID2=[9:11,14:22]; %%%%%%%% rotate chi1
    case 'Leu'
        numAtom =25;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 23];
        iChi1Array = [4,5,8,9];
        moveAtomID2=[9:11,14:22];
    case 'Val'
        numAtom = 22;
        Phi_ind = [2 4 5 6];
        Psi_ind = [ 4 5 6 20];
        iChi1Array = [4,5,8,9];
        moveAtomID2 = [9,10,13,14,15,16,17,18,19];
    case 'Phe'
        numAtom = 26;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 24];
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,10];
        moveAtomID=[10:14,19:23]; %%%%%%%% rotate chi2
        moveAtomID2=[9:14,17:23]; %%%%%%%% rotate chi1
        
    case 'Trp'
        numAtom = 30;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 28];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:17,20:27]; %%%%%%%% rotate chi1
    case 'Tyr'
        numAtom = 27;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 25];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:15,18:24]; %%%%%%%% rotate chi1
        
    case 'Asn'
        numAtom = 20;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 18];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:11,14:17]; %%%%%%%% rotate chi1
        
    case 'Cys'
        numAtom = 17;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 15];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9,12:14]; %%%%%%%% rotate chi1
    case 'Glu'
        numAtom =21;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 19];
        iChi1Array = [4,5,8,9];
        moveAtomID2=[9:12,15:18];
    case 'Met'
        numAtom = 23;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8,9,10,11];
        HG_Array_1 = [9,10,11,18];
        moveAtomID = [10,11,16,17,18,19,20];
        moveAtomID2 = [9,10,11,14,15,16,17,18,19,20];
        moveAtomID3 = [11,18,19,20];
        moveAtomID_HG1 = [18:20];
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 21];
        numAtom = 23;
    case 'Ser'
        numAtom =17;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 15];
        iChi1Array = [4,5,8,9];
        moveAtomID2=[9,12:14];
    case 'Thr'
        numAtom = 20;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 18];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:10,13:17]; %%%%%%%% rotate chi1
        
    case 'Asp'
        numAtom = 18;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 16];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:11,14:15]; %%%%%%%% rotate chi1
    case 'Gln'
        numAtom = 23;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 21];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:12,15:20]; %%%%%%%% rotate chi1
    case 'Arg'
        numAtom = 30;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 28];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:14,17:27]; %%%%%%%% rotate chi1
    case 'His'
        numAtom = 23;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 21];
        iChi1Array = [4,5,8,9];
        moveAtomID2=[9:13,16:20];
    case 'Lys'
        numAtom = 28;
        Phi_ind = [2 4 5 6];
        Psi_ind = [4 5 6 26];
        iChi1Array=[4,5,8,9];
        moveAtomID2=[9:12,15:25]; %%%%%%%% rotate chi1
        
    case 'Gly'
        fprintf('Not yet supported\n' );
    case 'Pro'
        fprintf('Not yet supported\n' );
    otherwise
        fprintf('Invalid amino acid\n' );
end
if next_pro == 1
    numAtom = numAtom-1;
end
%% Save dipeptide coordinates
XtalPosition = cell2mat(allDipeptide(:,8:10));
Position = XtalPosition;
Position_1000(:,1:3) = Position;

%Calculate Phi/Psi of original dipeptide
% And convert them to -180:180
Phi=calcDA2(Phi_ind,Position);
Psi=calcDA2(Psi_ind,Position);
if Phi > 180
    Phi = Phi - 360;
end
if Psi > 180
    Psi = Psi - 360;
end


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

i = 1;
for multi_pdb = 2:1000 %Make dipeptides with new bl/ba
    
    %Grab new coordinates
    new_res = Resi_coordinates(:,[(i-1)*3+1 : i*3]);
    if next_pro ==1
        new_res = new_res(1:size(new_res,1)-1, :);
    end
    
    %Move Ca to orgin
    TempPosition=Position-repmat(Position(5,:),numAtom,1);

    % Calculate what dihedral C-N-Ca-Cb angle should be
    InitAng=calcDA2([2 4 5 8],TempPosition);
    InitAng=mod(real(InitAng),360);
    
    %Calculate Chi1
    InitChi1=calcDA2(iChi1Array,TempPosition);
    InitChi1=mod(real(InitChi1),360);

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
    
    %% Orient original coordinates
    %Rotate everything around y axis until N is on x-y plane
    rot_Ny =  atan2d(TempPosition(4,3), TempPosition(4,1));
    R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
    TempPosition = (R*TempPosition')';
    
    
    %Rotate everything around z axis until N is on x axis
    rot_Ny =  -atan2d(TempPosition(4,2), TempPosition(4,1));
    R =    [ cosd(rot_Ny) ,-sind(rot_Ny),0;sind(rot_Ny), cosd(rot_Ny),0;0 0 1];
    TempPosition = (R*TempPosition')';
    if TempPosition(4,1) < 0
        rot_Ny =  180;
        R =    [ cosd(rot_Ny) ,0, sind(rot_Ny); 0 1 0; -sind(rot_Ny) 0 cosd(rot_Ny)];
        TempPosition = (R*TempPosition')';
    end
    
    %Now N is on the x axis, move Cb to the x,y plane, rotate on x
    rot_Ny =  -atan2d(TempPosition(8,3), TempPosition(8,2));
    R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
    TempPosition = (R*TempPosition')';
    if TempPosition(8,2) < 0
        rot_Ny = 180;
        R =    [ 1 0 0; 0 cosd(rot_Ny) -sind(rot_Ny); 0 sind(rot_Ny) cosd(rot_Ny)];
        TempPosition = (R*TempPosition')';
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