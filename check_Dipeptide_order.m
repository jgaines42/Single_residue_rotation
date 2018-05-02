%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(residue, Resi_name)
%
% Checks the order of dipeptide atoms and returns them in a standard order
%
% Input: 
%   residue: The cell array of the residue (not the Dipeptide)
%   Resi_name: Residue name (3 letter code)
%
% Output: 
%   new_residue: The residue in the correct order
%   correct_now: 1 if now in the correct order, 0 if not
%   ind: indexing of the new order (from the original order)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [new_residue, correct_now, ind] = check_Dipeptide_order(residue, Resi_name)
new_residue = residue;
correct_now = 0;
ind = 0;
switch (Resi_name)
    case 'Cys'
        if size(residue,1) == 17-6
            order = { 'N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', '1HB', '2HB', 'HG'};
          
            [y, ind] = ismember(order, residue(:,2));
           
            if y
                correct_now = 1;
                 new_residue = residue(ind,:);
            else
                order = {'N';'CA';'C';'O';'CB';'SG';'H';'HA';'HB2';'HB3';'HG'};
                [y,ind]= ismember(order, residue(:,2));
                if y
                    new_residue = residue(ind,:);
                    correct_now = 1;
                end
            end
        elseif size(residue,1) == 16-6
            order = {'N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', '1HB', '2HB'};
            [y, ind] = ismember(order, residue(:,2));
           
            if y
                 new_residue = residue(ind,:);
                correct_now = 1;
            end
        end
        
    case 'Ile'
        order = { 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', '1HG1','2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order = {'N';'CA';'C';'O';'CB';'CG1';'CG2';'CD1';'H';'HA';'HB';'HG12';'HG13';'HG21';'HG22';'HG23';'HD11';'HD12';'HD13'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
       
        end
        
    case 'Leu'
         order = { 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', '1HB', '2HB', 'HG', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
          else
            order =   {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'H';'HA';'HB2';'HB3';'HG';'HD11';'HD12';'HD13';'HD21';'HD22';'HD23'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
       
        end   
            
          
       
        
    case 'Met'
         order = { 'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', '1HB', '2HB','1HG', '2HG', '1HE', '2HE', '3HE'};
        
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order =   {'N';'CA';'C';'O';'CB';'CG';'SD';'CE';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HE1';'HE2';'HE3'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            else
                order =   {'N';'CA';'C';'O';'CB';'CG';'SD';'CE';'H';'HA';'HB1';'HB2';'HG1';'HG2';'HE1';'HE2';'HE3'};
                [y,ind]= ismember(order, residue(:,2));
                if y
                    new_residue = residue(ind,:);
                    correct_now = 1;
                end

            end
        end
        
    case 'Phe'
        order = { 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1','CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', '1HB', '2HB','HD1', 'HD2', 'HE1', 'HE2', 'HZ'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order =     {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'CE1';'CE2';'CZ';'H';'HA';'HB2';'HB3';'HD1';'HD2';'HE1';'HE2';'HZ'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
        end
        
        
       
    case 'Ser'
        order = { 'N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', '1HB', '2HB','HG'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order =    {'N';'CA';'C';'O';'CB';'OG';'H';'HA';'HB2';'HB3';'HG'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
        end
        
        
      
    case 'Thr'

        order = { 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'H', 'HA', 'HB', 'HG1', '1HG2', '2HG2', '3HG2'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order =   {'N';'CA';'C';'O';'CB';'OG1';'CG2';'H';'HA';'HB';'HG1';'HG21';'HG22';'HG23'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
        end
        
    case 'Tyr'
        
         order = { 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1','CD2', 'CE1', 'CE2', 'CZ', 'OH','H', 'HA', '1HB', '2HB','HD1', 'HD2', 'HE1', 'HE2', 'HH'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order =   {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'CE1';'CE2';'CZ';'OH';'H';'HA';'HB2';'HB3';'HD1';'HD2';'HE1';'HE2';'HH'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
            
        end
    case 'Val'
         order = { 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', '1HG1','2HG1','3HG1', '1HG2', '2HG2', '3HG2'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        else
            order = {'N';'CA';'C';'O';'CB';'CG1';'CG2';'H';'HA';'HB';'HG11';'HG12';'HG13';'HG21';'HG22';'HG23'};
            [y,ind]= ismember(order, residue(:,2));
            if y
               new_residue = residue(ind,:);
                correct_now = 1;
            end
        end
    case 'Ala'
         order = {'N';'CA';'C';'O';'CB';'H';'HA';'HB1';'HB2';'HB3'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
%         else
%             order = {'N';'CA';'C';'O';'CB';'CG1';'CG2';'H';'HA';'HB';'HG11';'HG12';'HG13';'HG21';'HG22';'HG23'};
%             [y,ind]= ismember(order, residue(:,2));
%             if y
%                new_Dipeptide = residue(ind,:);
%                 correct_now = 1;
%             end
        end
    case 'Trp'
        order = {'N';'CA';'C';'O';'CB';'CG';'CD1';'CD2';'NE1';'CE2';'CE3';'CZ2';'CZ3';'CH2';'H';'HA';'HB2';'HB3';'HD1';'HE1';'HE3';'HZ2';'HZ3';'HH2'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        end
    
    case 'His'
        order =  {'N';'CA';'C';'O';'CB';'CG';'ND1';'CD2';'CE1';'NE2';'H';'HA';'HB2';'HB3';'HD2';'HE1'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        end
    case 'Lys'
        order =  {'N';'CA';'C';'O';'CB';'CG';'CD';'CE';'NZ';'H';'HA';'HB2';'HB3';'HG2';'HG3';'HD2';'HD3';'HE2';'HE3';'HZ1';'HZ2';'HZ3'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        end
    case 'Asp'
        order = {'N';'CA';'C';'O';'CB';'CG';'OD1';'OD2';'H';'HA';'HB2';'HB3'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        end
    case 'Glu'
        order = {'N';'CA';'C';'O';'CB';'CG';'CD';'OE1';'OE2';'H';'HA';'HB2';'HB3';'HG2';'HG3'};
        [y, ind] = ismember(order, residue(:,2));
        if y
            new_residue = residue(ind,:);
            correct_now = 1;
        end
    otherwise
        ind = 0;
        correct_now = 0;
    
    
end

end