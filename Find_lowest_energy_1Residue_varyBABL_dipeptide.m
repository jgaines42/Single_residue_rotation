%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [min_energy,dihedrals] = Find_lowest_energy_1Residue_varyBABL(Position, next_pro, energy_cutoff, Amino_Acid, Init_dihedrals)
%
% Rotates a residue and finds the lowest energy state. Returns a list of all
% dihedral angles with that energy
%
% Input:
%   Position: Position coordinates for the dipeptide
%   res_name: Residue name (3 letter)
%   Atom_sizes: Atom sizes of the dipeptide atoms
%   next_pro: 1 if next residue is a proline
%   energy_cutoff: cutoff for energy search, finds minimum energy below
%   this value
%   Amino_Acid: Structure of data related to the amino acid
%   Init_dihedrals: Structure of initial dihedral angle values of the amino
%   acid
%
% Output:
%   min_energy: Lowest energy of the residue
%   dihedrals: All dihedral angles with that energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [min_energy,dihedrals] = Find_lowest_energy_1Residue_varyBABL_dipeptide(Position, energy_cutoff, clashLists, Amino_Acid, Init_dihedrals)

XtalPos = Position;

min_energy = energy_cutoff;
dihedrals = [];

%Extract data from structures (cause its faster this way)
DOF = Amino_Acid.DOF;
CH3 = Amino_Acid.CH3;
OH = Amino_Acid.OH;

iChi1Array = Amino_Acid.iChi1Array;
moveAtomID2 = Amino_Acid.moveAtomID2;
c1_Clash = clashLists.c1_Clash;
InitChi1 = Init_dihedrals.InitChi1;

if DOF > 1
    iChi2Array = Amino_Acid.iChi2Array;
    moveAtomID = Amino_Acid.moveAtomID;
    c2_Clash = clashLists.c2_Clash;
    InitChi2 = Init_dihedrals.InitChi2;
    if DOF > 2
        iChi3Array = Amino_Acid.iChi3Array;
        moveAtomID3 = Amino_Acid.moveAtomID3;
        c3_Clash = clashLists.c3_Clash;
        InitChi3 = Init_dihedrals.InitChi3;
    end
end
if CH3 >= 1
    HG_Array_1 = Amino_Acid.HG_Array_1;
    moveAtomID_HG1 = Amino_Acid.moveAtomID_HG1;
    HG1_Clash = clashLists.HG1_Clash;
    InitChi_HG1 = Init_dihedrals.InitChi_HG1;
    
    if CH3 == 2
        HG_Array_2 = Amino_Acid.HG_Array_2;
        moveAtomID_HG2 = Amino_Acid.moveAtomID_HG2;
        HG2_Clash = clashLists.HG2_Clash;
        InitChi_HG2 = Init_dihedrals.InitChi_HG2;
        
    end
end
if OH == 1
    iOHArray= Amino_Acid.iOHArray;
    moveAtomOH = Amino_Acid.moveAtomOH;
    OH_Clash = clashLists.OH_Clash;
    InitOH = Init_dihedrals.InitOH;
    
end

numAtom = size(Position,1);

%% Start rotations
subtract_array_1 = repmat(Position(iChi1Array(2),:),numAtom,1);
delta_term_1 =  pi*sign(InitChi1)*InitChi1/180;

for chi1 = 1:72 % Loop through chi1 sampling every 5 degrees
    Position=XtalPos;
    setChi1 = chi1*5;
    Position = Rotate_DA(Position, setChi1, subtract_array_1, delta_term_1, iChi1Array, moveAtomID2);
    chi1_energy = get_energy_wProtein(c1_Clash, Position,  [], []);
    
    %% Only enter chi2 loop if the energy is lower than the minimum energy so far
    if DOF >=2 && chi1_energy <= min_energy
        Pos_b4_Chi2 = Position;
        subtract_array_2 = repmat(Position(iChi2Array(2),:),numAtom,1);
        delta_term_2 =  pi*sign(InitChi2)*InitChi2/180;
        
        for chi2 = 1:Init_dihedrals.max_Chi2
            Position=Pos_b4_Chi2;
            setChi2 = chi2*5;
            Position = Rotate_DA(Position, setChi2, subtract_array_2, delta_term_2, iChi2Array, moveAtomID);
            [chi2_energy] = get_energy_wProtein(c2_Clash, Position,  [], []);
            
            total_energy = chi1_energy + chi2_energy;
            
            %% Only enter chi3 loop if the energy is lower than the minimum energy so far
            if DOF >= 3 && total_energy <= min_energy
                Pos_b4_Chi3 = Position;
                subtract_array_3 = repmat(Position(iChi3Array(2),:),numAtom,1);
                delta_term_3 =  pi*sign(InitChi3)*InitChi3/180;
                
                for chi3 = 1:72
                    Position=Pos_b4_Chi3;
                    setChi3 = chi3*5;
                    Position = Rotate_DA(Position, setChi3, subtract_array_3, delta_term_3, iChi3Array, moveAtomID3);
                    [chi3_energy] = get_energy_wProtein(c3_Clash, Position,  [], []);
                    total_energy = chi1_energy + chi2_energy + chi3_energy;
                    
                    %% Check for methyl group at end and rotate
                    if DOF == 3 && total_energy <= min_energy
                        this_dihedral = [setChi1, setChi2, setChi3];
                        [dihedrals, min_energy] = check_end_groups(dihedrals, min_energy,total_energy, this_dihedral, Amino_Acid,Init_dihedrals, Position, clashLists);
                        
                    end
                end%Chi3 loop
            end
            
            %% Check for methyl group at end and rotate
            if DOF == 2 && total_energy <= min_energy
                
                this_dihedral = [setChi1, setChi2];
                [dihedrals, min_energy] = check_end_groups(dihedrals, min_energy,total_energy, this_dihedral, Amino_Acid,Init_dihedrals, Position, clashLists);
                
            end
            
            
        end
        
    end
    
    %% Check for methyl group at end and rotate
    if DOF == 1 && chi1_energy <= min_energy
        total_energy = chi1_energy;
        this_dihedral = [setChi1];
        [dihedrals, min_energy] = check_end_groups(dihedrals, min_energy,total_energy, this_dihedral, Amino_Acid,Init_dihedrals, Position, clashLists);
        
        
    end
end

end
