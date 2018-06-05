%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function tempModel2 = add_sizes_protein(tempModel2,size_number)
%
% Adds radii to all atoms in the protein, stores in column 12
% 
% Input:
% tempModel2: The protein structure
% size_number: Which interation of atom sizes to use
%
% Returns:
% tempModel2 with sizes stored in column 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function tempModel2 = add_sizes_protein(tempModel2,size_number)

 tempModel2(:,12) = {0};
 
% Order of atoms:
% 1     2   3  4    5       6       7   8  9      10    11    12      13     14         15    16 
% H_C, H_N, N, O, C(sp3), C(sp2), O(H), S, H_S, C(Ar), H(Ar), O(Ar),  C(NO), O(CNO),  C(OH), H_O 
if size_number == 6
    radii_array = [1.1 1.1 1.3 1.4 1.5 1.3 1.4 1.75 1.1 1.50 1.1 1.4 1.3 1.4 1.5 1.1 1.3 1.3];
elseif size_number == 7
    radii_array = [1.1 1.0 1.3 1.4 1.5 1.3 1.4 1.75 1.1 1.50 1.1 1.4 1.3 1.4 1.5 1.0 1.3 1.3];
elseif size_number == 9
    radii_array = [1.1 1.0 1.3 1.4 1.5 1.3 1.4 1.75 1.0 1.45 1.05 1.4 1.3 1.4 1.5 1.0 1.3 1.3];
elseif size_number == 71
    radii_array = [1.1 1.1 1.3 1.4 1.5 1.3 1.4 1.75 1.1 1.45 1.05 1.4 1.3 1.4 1.5 1.0 1.3 1.3];
elseif size_number == 72
    radii_array = [1.1 1.0 1.3 1.4 1.5 1.3 1.4 1.75 1.1 1.45 1.1 ];
elseif size_number == 73
    radii_array = [1.1 1.1 1.3 1.4 1.5 1.3 1.4 1.75 1.0 1.45 1.05 ];
elseif size_number == 74
    radii_array = [1.1 1.1 1.3 1.4 1.5 1.3 1.4 1.75 1.0 1.45 1.05 0 0 0  1.5 1.0];
elseif size_number == 1
    radii_array = [0    0    1.64    1.42    1.88  1.61 1.46 1.77 0 1.45 1.05 1.4 0 0 1.425];
end

Ile_atoms = {'N';'CA';'C';'O';'CB';'CG1';'CG2';'CD1';'H';'HA';'HB';'HG12';'HG13';'HG21';'HG22';'HG23';'HD11';'HD12';'HD13'}';
Leu_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'H';'HA';'HB2';'HB3';'HG';'HD11';'HD12';'HD13';'HD21';'HD22';'HD23'}';
Phe_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'CE1';'CE2';'CZ';'H';'HA';'HB2';'HB3';'HD1';'HD2';'HE1';'HE2';'HZ' }';
Val_atoms = {'N';'CA';'C';'O';'CB';'CG1';'CG2';'H';'HA';'HB';'HG11';'HG12';'HG13';'HG21';'HG22';'HG23' }';
His_atoms = {'N';'CA';'C';'O';'CB';'CG';'ND1';'CD2';'CE1';'NE2';'H';'HA';'HB2';'HB3';'HD2';'HE1'; 'HD1'; 'HE2'}';
Gln_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD';'OE1';'NE2';'H';'HA'; 'HB2';'HB3'; 'HG2'; 'HG3'; 'HE21'; 'HE22'}';
Arg_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD';'NE';'CZ';'NH1';'NH2';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HD2';'HD3';'HE';'HH11';'HH12';'HH21';'HH22'}';
Trp_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'NE1';'CE2';'CE3';'CZ2';'CZ3';'CH2';'H';'HA'; 'HB2'; 'HB3';'HD1';'HE1';'HE3';'HZ2';'HZ3';'HH2';}';
Asp_atoms = {'N';'CA';'C';'O';'CB';'CG';'OD1';'OD2';'H';'HA'; 'HB2'; 'HB3'}';
Asn_atoms = {'N';'CA';'C';'O';'CB';'CG';'OD1';'ND2';'H';'HA';'HB2';'HB3';'HD21';'HD22'}';
Cys_atoms = {'N';'CA';'C';'O';'CB';'SG';'H';'HA';'HB2';'HB3';'HG';}';
Met_atoms = {'N';'CA';'C';'O';'CB';'CG';'SD';'CE';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HE1';'HE2';'HE3'}';
Ser_atoms = {'N';'CA';'C';'O';'CB';'OG';'H';'HA';'HB2'; 'HB3';'HG'}';
Thr_atoms = {'N';'CA';'C';'O';'CB';'OG1';'CG2';'H';'HA';'HB';'HG1';'HG21';'HG22';'HG23'}';
Tyr_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'CE1';'CE2';'CZ';'OH';'H';'HA';'HB2'; 'HB3';'HD1';'HD2';'HE1';'HE2';'HH'}';
Glu_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD';'OE1';'OE2';'H';'HA';'HB2';'HB3';'HG2';'HG3'}';
Lys_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD';'CE';'NZ';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HD2';'HD3';'HE2';'HE3';'HZ1';'HZ2';'HZ3'}';
Pro_atoms = {'N';'CA';'C';'O';'CB';'CG';'CD';'HA';'HB2';'HB3';'HG2';'HG3';'HD2';'HD3'}';
Ala_atoms = {'N';'CA';'C';'O';'CB';'H';'HA';'HB1';'HB2';'HB3';}';
Gly_atoms = {'N';'CA';'C';'O';'H'; 'HA2'; 'HA3'}';
Mse_atoms = {'N';'CA';'C';'O';'CB';'CG';'SE';'CE';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HE1';'HE2';'HE3'}';
all_res_names = {'Ala', 'Ile', 'Leu', 'Phe', 'Val', 'His', 'Gln', 'Arg', 'Trp', 'Asp', 'Asn', 'Cys', 'Met', 'Ser', 'Thr', 'Tyr', 'Glu', 'Lys', 'Pro', 'Gly', 'Mse'};
for amino_acids = 1:20
    res_name = all_res_names{amino_acids};
    ind0 = strcmpi(tempModel2(:,4), res_name);
    all_this_res = tempModel2(ind0,:);
    
    switch(res_name)
        case 'Ala'
            size_index = [5;6;4;3;5;6;4;5;2;1;1;1;1;3;5;2];
            res_atoms = Ala_atoms;
        case 'Ile'
            size_index = [5;6;4;3;5;6;4;5;5;5;5;2;1;1;1;1;1;1;1;1;1;1;3;5;2];
            res_atoms = Ile_atoms;
        case 'Leu'
            size_index = [5;6;4;3;5;6;4;5;5;5;5;2;1;1;1;1;1;1;1;1;1;1;3;5;2];
            res_atoms = Leu_atoms;
        case 'Val'
            size_index = [5;6;4;3;5;6;4;5;5;5;2;1;1;1;1;1;1;1;1;3;5;2];
            res_atoms =  Val_atoms;
        case 'Phe'
            res_atoms = Phe_atoms;
            size_index = [5;6;4;3;5;6;4;5;10;10;10;10;10;10;2;1;1;1;11;11;11;11;11;3;5;2];
        case 'Trp'
            res_atoms = Trp_atoms;
            size_index = [5;6;4;3;5;6;4;5;10;10;10;3;10;10;10;10;10;2;1;1;1;11;2;11;11;11;11;3;5;2];
        case 'Tyr'
            res_atoms = Tyr_atoms;
            size_index = [5; 6; 4 ;3; 5; 6 ;4 ;5; 10;10;10;10;10;10; 12; 2; 1; 1 ;1 ; 11;11;11;11; 2 ;3; 5 ;2];
        case 'Asn'
            res_atoms = Asn_atoms;
            size_index = [5; 6; 4 ;3; 5; 6 ;4; 5; 13; 14; 17; 2;1;1;1;2;2;3;5;2];
        case 'Cys'
            res_atoms = Cys_atoms;
            size_index = [5;6;4;3;5;6;4;5;8;2;1;1;1;9;3;5;2];
        case 'Glu'
            res_atoms = Glu_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;13;14;14;2;1;1;1;1;1;3;5;2];
        case 'Met'
            res_atoms = Met_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;8;5;2;1;1;1;1;1;1;1;1;3;5;2];
        case 'Ser'
            res_atoms = Ser_atoms;
            size_index = [5;6;4;3;5;6;4;15;7;2;1;1;1;16;3;5;2];
        case 'Thr'
            res_atoms = Thr_atoms;
            size_index = [5;6;4;3;5;6;4;15;7;5;2;1;1;16;1;1;1;3;5;2];
        case 'Asp'
            res_atoms = Asp_atoms;
            size_index = [5;6;4;3;5;6;4;5;13;14;14;2;1;1;1;3;5;2];
        case 'Gln'
            res_atoms = Gln_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;13;14;17;2;1;1;1;1;1;2;2;3;5;2];
        case 'Arg'
            res_atoms = Arg_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;5;3;18;17;17;2;1;1;1;1;1;1;1;2;2;2;2;2;3;5;2];
        case 'His'
            res_atoms = His_atoms;
            size_index = [5;6;4;3;5;6;4;5;10;3;10;10;3;2;1;1;1;2;2;3;5;2];
        case 'Lys'
            res_atoms = Lys_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;5;5;3;2;1;1;1;1;1;1;1;1;1;2;2;2;3;5;2];
        case 'Gly'
            res_atoms = Gly_atoms;
            size_index = [5;6;4;3;5;6;4;2;1;1;3;5;2];
        case 'Pro'
            res_atoms = Pro_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;5;1;1;1;1;1;1;1;3;5;2];
        case 'Mse'
            res_atoms = Met_atoms;
            size_index = [5;6;4;3;5;6;4;5;5;8;5;2;1;1;1;1;1;1;1;1;3;5;2];
        otherwise
            fprintf('Invalid amino acid\n' );
    end
    size_index = size_index(4:size(size_index,1)-3);
    if res_name == 'His'
        size_index =[ size_index; 1;1];
    end
    Atom_sizes =  radii_array(1, size_index)';
    
    for i = 1:size(res_atoms, 2)
        if i ==6
            here = 1;
        end
        ind0 = strcmp(all_this_res(:,2), res_atoms(i));
        tempModel2(cell2mat(all_this_res(ind0,1)),12) = num2cell(Atom_sizes(i));
    end
end
ind0 = ismember(tempModel2(:,2), {'H1', 'H2' ,'H3'});
tempModel2(ind0,12) = {1.1};

ind0 = ismember(tempModel2(:,2), {'OXT'});
tempModel2(ind0,12) = {1.3};
end