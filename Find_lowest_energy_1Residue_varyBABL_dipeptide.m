%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [min_energy,dihedrals] = Find_lowest_energy_1Residue_varyBABL(Position, res_name, Atom_sizes, rest_of_pro, next_pro, energy_cutoff)
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
%   this
%
% Output:
%   min_energy: Lowest energy of the residue
%   dihedrals: All dihedral angles with that energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [min_energy,dihedrals] = Find_lowest_energy_1Residue_varyBABL_dipeptide(Position, res_name, Atom_sizes, next_pro, energy_cutoff)

XtalPos = Position;

%energy_cutoff = 1;
min_energy = energy_cutoff;
dihedrals = [];

%Run setup script to get all of the variables for the particular amino acid
switch_residue_setup;

set_up_clashArrays_dipeptide;

numAtom = size(Position,1);

subtract_array_1 = repmat(Position(iChi1Array(2),:),numAtom,1);
delta_term_1 =  pi*sign(InitChi1)*InitChi1/180;

for chi1 = 1:72
    Position=XtalPos;
    setChi1 = chi1*5;
    Position = Rotate_DA(Position, setChi1, subtract_array_1, delta_term_1, iChi1Array, moveAtomID2);
    [chi1_energy] = get_energy_wProtein(c1_Clash, Position,  [], []);
    
    if DOF >=2 && chi1_energy <= min_energy
        Pos_b4_Chi2 = Position;
        subtract_array_2 = repmat(Position(iChi2Array(2),:),numAtom,1);
        delta_term_2 =  pi*sign(InitChi2)*InitChi2/180;
        for chi2 = 1:max_Chi2
            Position=Pos_b4_Chi2;
            setChi2 = chi2*5;
            Position = Rotate_DA(Position, setChi2, subtract_array_2, delta_term_2, iChi2Array, moveAtomID);
            [chi2_energy] = get_energy_wProtein(c2_Clash, Position,  [], []);
            
            total_energy = chi1_energy + chi2_energy;
            
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
                    
                    if DOF == 3 && total_energy <= min_energy
                        if CH3 == 1
                            subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                            delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                            min_energy_CH3  = rotate_CH3_group_1_find_lowest(Position, HG1_Clash, [], [], subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
                            
                            total_energy = total_energy + min_energy_CH3;
                        end
                        if total_energy < min_energy
                            min_energy = total_energy;
                            dihedrals = [setChi1, setChi2, setChi3];
                        elseif total_energy == min_energy
                            dihedrals = [dihedrals;setChi1, setChi2, setChi3];
                        end
                        
                    end
                end%Chi3 loop
            end
            
            if DOF == 2 && total_energy <= min_energy
                if CH3 == 1
                    
                    subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                    delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                    min_energy_CH3  = rotate_CH3_group_1_find_lowest(Position, HG1_Clash, [], [], subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
                    total_energy = total_energy + min_energy_CH3;
                    
                elseif CH3 == 2
                    
                    subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                    delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                    delta_term_HG2 =  pi*sign(InitChi_HG2)*InitChi_HG2/180;

                    min_energy_CH3 = rotate_CH3_group_2_find_lowest(Position, HG1_Clash, [], HG2_Clash, [], [], subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1, HG_Array_2, moveAtomID_HG2, numAtom,delta_term_HG2,min_energy, total_energy);
                    
                    total_energy = total_energy + min_energy_CH3;
                end
                if total_energy < min_energy
                    min_energy = total_energy;
                    dihedrals = [setChi1, setChi2];
                elseif total_energy == min_energy
                    dihedrals = [dihedrals;setChi1, setChi2];
                end
            end
            
            
        end
        
    end
    
    if DOF == 1 && chi1_energy <= min_energy
        total_energy = chi1_energy;
        if CH3 == 1 &&  OH == 0 && DOF == 1
            
            subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
            delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
            min_energy_CH3  = rotate_CH3_group_1_find_lowest(Position, HG1_Clash, [], [], subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
            
            total_energy = chi1_energy + min_energy_CH3;
            
        elseif OH == 0 && DOF == 1 && CH3 == 2
            subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
            delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
            delta_term_HG2 =  pi*sign(InitChi_HG2)*InitChi_HG2/180;

            min_energy_CH3 = rotate_CH3_group_2_find_lowest(Position, HG1_Clash, [], HG2_Clash, [], [], subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1, HG_Array_2, moveAtomID_HG2, numAtom,delta_term_HG2,min_energy, chi1_energy);
            
            total_energy = chi1_energy + min_energy_CH3;
            
            %%%%%%%
        elseif DOF == 1 && OH == 1 && CH3 == 0
            Pos_b4_OH = Position;
            subtract_array_OH = repmat(Position(iOHArray(2),:),numAtom,1);
            delta_term_OH =  pi*sign(InitOH)*InitOH/180;
            min_energy_OH  = rotate_CH3_group_1_find_lowest(Position, OH_Clash, [], [], subtract_array_OH, delta_term_OH, iOHArray, moveAtomOH);
            total_energy = chi1_energy +min_energy_OH;
            
        elseif   DOF == 1 && OH == 1 && CH3 == 1
            Pos_b4_HG1 = Position;
            subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
            delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
            
            
            [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  [], []);
            [OH_energy] = get_energy_wProtein(OH_Clash, Position,  [], []);
            
            min_HG1 = HG1_energy;
            min_OH = OH_energy;
            min_OH_HG1 = min_HG1 + min_OH;
            
            HG1_val = 0;
            
            while min_OH_HG1 > 0 && HG1_val < 72
                HG1_val = HG1_val +1;
                setHG1 = HG1_val*5;
                Position=Pos_b4_HG1;
                Position = Rotate_DA(Position, setHG1, subtract_array_HG1, delta_term_HG1, HG_Array_1, moveAtomID_HG1);
                [HG1_energy] = get_energy_wProtein(HG1_Clash, Position,  [], []);
                
                
                if HG1_energy <= min_OH_HG1
                    Pos_b4_OH = Position;
                    subtract_array_OH = repmat(Position(iOHArray(2),:),numAtom,1);
                    delta_term_OH =  pi*sign(InitOH)*InitOH/180;
                    
                    [OH_energy] = get_energy_wProtein(OH_Clash, Position,  [], []);
                    
                    min_OH = OH_energy;
                    
                    OH_val = 0;
                    while min_OH>0 && OH_val < 72
                        OH_val = OH_val +1;
                        setOH = OH_val*5;
                        Position = Pos_b4_OH;
                        Position = Rotate_DA(Position, setOH, subtract_array_OH, delta_term_OH, iOHArray, moveAtomOH);
                        [temp_energy] = get_energy_wProtein(OH_Clash, Position,  [], []);
                        if temp_energy > 0 && temp_energy < min_OH
                            OH_energy = temp_energy;
                            min_OH = temp_energy;
                        end
                        if temp_energy == 0
                            OH_energy = 0;
                            min_OH = 0;
                        end
                    end
                    if min_OH + HG1_energy < min_OH_HG1
                        min_OH_HG1 = min_OH + HG1_energy;
                    end
                end
            end
            total_energy = chi1_energy + min_OH_HG1;
            
            
            %%%%%
        end
        if total_energy < min_energy
            min_energy = total_energy;
            dihedrals = [setChi1];
        elseif total_energy == min_energy;
            dihedrals = [dihedrals;setChi1];
        end
        
    end
end

end