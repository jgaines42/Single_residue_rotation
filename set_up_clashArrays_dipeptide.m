%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [Init_dihedrals, Clash_lists] = set_up_clashArrays_dipeptide( res_name,next_pro,Atom_sizes, Amino_acid, Position)
%
% Sets up the clash lists for each dihedral angle movement in the dipeptide
% and calculates the initial dihedral angles
%
% Input:
% res_name: 3 letter abbreviation of the amino acid
% next_pro: 1 if the next amino acid is Pro, 0 if it isn't
% Atom_sizes: The atoms sizes for the entire dipeptide
% Amino_acid: Amino_acid structure created by switch_residue_setup()
% Position: Positions of the atoms in the dipeptide
%
% Output:
% Init_dihedrals: Structure containing initial dihedral angle values
% Clash_lists: Structures containing all clash lists
%
% Notes:
% For each clashlist, only the atoms that move in that dihedral roation but
% not the next rotation are included. For example, if rotating Chi1 moves
% atoms 9 and 10 while rotation Chi2 moves atom 10, clashes with atom 9
% will be included in c1_clashList but not atom 10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Init_dihedrals, Clash_lists] = set_up_clashArrays_dipeptide( res_name,next_pro,Atom_sizes, Amino_acid, Position)
CH3 = Amino_acid.CH3;
OH = Amino_acid.OH;
DOF = Amino_acid.DOF;
Init_dihedrals.max_Chi2 = 0;
if DOF >= 3
    load(strcat('clash_folder/chi3_', res_name, '_clash.mat'));
    c3_clashList = n;
    if next_pro == 1
        ind0 = c3_clashList(:,2) == numAtom;
        c3_clashList = c3_clashList(~ind0,:);
    end
    
     Clash_lists.c3_Clash = [c3_clashList, (Atom_sizes(c3_clashList(:,1),:)+Atom_sizes(c3_clashList(:,2),:))];
    InitChi3=calcDA2( Amino_acid.iChi3Array,Position); %%%%%%%%
    Init_dihedrals.InitChi3=mod(real(InitChi3),360);
    Clash_lists.c3_Clash(:,3)=Clash_lists.c3_Clash(:,3).^2;
end

if DOF >= 2
    load(strcat('clash_folder/chi2_', res_name, '_clash.mat'));
    c2_clashList = n;
    if next_pro == 1
        ind0 = (c2_clashList(:,2) == numAtom);
        c2_clashList = c2_clashList(~ind0,:);
    end
     Clash_lists.c2_Clash = [c2_clashList, (Atom_sizes(c2_clashList(:,1),:)+Atom_sizes(c2_clashList(:,2),:))];
    Clash_lists.c2_Clash(:,3)= Clash_lists.c2_Clash(:,3).^2;
    
    
    InitChi2=calcDA2( Amino_acid.iChi2Array,Position); %%%%%%%%
   Init_dihedrals.InitChi2=mod(real(InitChi2),360);
     
   % Because Phe, tyr and Asp are symetric around chi2, change the max to
   % sample
   Init_dihedrals.max_Chi2 = 72;
    if strcmp(res_name,'Phe') || strcmp(res_name, 'Tyr') || strcmp(res_name, 'Asp')
        Init_dihedrals.max_Chi2 = 72/2;
    end
end
if DOF >= 1
    load(strcat('clash_folder/chi1_', res_name, '_clash.mat'));
    c1_clashList = n;
    if next_pro == 1
        ind0 = c1_clashList(:,2) == numAtom;
        c1_clashList = c1_clashList(~ind0,:);
    end
     Clash_lists.c1_Clash = [c1_clashList, (Atom_sizes(c1_clashList(:,1),:)+Atom_sizes(c1_clashList(:,2),:))];
    
    
    InitChi1=calcDA2( Amino_acid.iChi1Array,Position);
    Init_dihedrals.InitChi1=mod(real(InitChi1),360);
  
    Clash_lists.c1_Clash(:,3) =  Clash_lists.c1_Clash(:,3).^2;
    
end
if CH3 == 2
    load(strcat('clash_folder/HD2_', res_name, '_clash.mat'));
    HG2_clashList = n;
    if next_pro == 1
        ind0 = HG2_clashList(:,2) == numAtom;
        HG2_clashList = HG2_clashList(~ind0,:);
    end
    Clash_lists.HG2_Clash = [HG2_clashList, (Atom_sizes(HG2_clashList(:,1),:)+Atom_sizes(HG2_clashList(:,2),:))];
    Clash_lists.HG2_Clash(:,3) =  Clash_lists.HG2_Clash(:,3).^2;
    
    InitChi_HG2=calcDA2( Amino_acid.HG_Array_2,Position);
    Init_dihedrals.InitChi_HG2=mod(real(InitChi_HG2),360);
    
end
if CH3 >=1
    load(strcat('clash_folder/HD1_', res_name, '_clash.mat'));
    HG1_clashList = n;
    if next_pro == 1
        ind0 = HG1_clashList(:,2) == numAtom;
        HG1_clashList = HG1_clashList(~ind0,:);
    end
    Clash_lists.HG1_Clash = [HG1_clashList, (Atom_sizes(HG1_clashList(:,1),:)+Atom_sizes(HG1_clashList(:,2),:))];
    Clash_lists.HG1_Clash(:,3) =  Clash_lists.HG1_Clash(:,3).^2;
 
    InitChi_HG1=calcDA2( Amino_acid.HG_Array_1,Position); 
    Init_dihedrals.InitChi_HG1=mod(real(InitChi_HG1),360);
    
end
if OH == 1
    if res_name == 'Cys'
        load(strcat('clash_folder/SH_', res_name, '_clash.mat'));
    else
        load(strcat('clash_folder/OH_', res_name, '_clash.mat'));
    end
    OH_clashList = n;
    if next_pro == 1
        ind0 = OH_clashList(:,2) == numAtom;
        OH_clashList = OH_clashList(~ind0,:);
    end
    Clash_lists.OH_Clash = [OH_clashList, (Atom_sizes(OH_clashList(:,1),:)+Atom_sizes(OH_clashList(:,2),:))];
    Clash_lists.OH_Clash(:,3) =  Clash_lists.OH_Clash(:,3).^2;
    
    InitOH=calcDA2( Amino_acid.iOHArray,Position); %InitChi1
    Init_dihedrals.InitOH=mod(real(InitOH),360); 
end

Clash_lists.protein_clash_c1 = [];
Clash_lists.protein_clash_c2 = [];
Clash_lists.protein_clash_c3 = [];
Clash_lists.protein_clash_HG1 = [];
Clash_lists.protein_clash_HG2 = [];
Clash_lists.protein_clash_OH = [];

end
