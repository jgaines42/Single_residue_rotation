function [] = process_single_cluster(what_res, which_envir, which_database)

if strcmp(which_envir, 'pro')
    sub_name = '';
else
    sub_name ='_dipeptide';
end

Resi_name = upper(what_res);
if strcmp(which_database,'Dun')
    load(strcat('/ysm-gpfs/home/jcg72/Rotation_SASA/Dun_1A/', what_res, '_locations.mat'));
    data_folder = strcat('/ysm-gpfs/home/jcg72/Rotation_SASA/Dun_',what_res,'/');

else
    load(strcat('/ysm-gpfs/home/jcg72/Rotation_SASA/PPI/', what_res, '_locations.mat'));
    data_folder = strcat('/ysm-gpfs/home/jcg72/Rotation_SASA/PPI/PPI_',what_res,'/');

end
switch what_res
    case 'Ile'
        all_loc = all_Ile;
    case 'Leu'
        all_loc = all_Leu;
    case 'Met'
        all_loc = all_Met;
    case 'Phe'
        all_loc = all_Phe;
    case 'Ser'
        all_loc = all_Ser;
    case 'Thr'
        all_loc = all_Thr;
    case 'Trp'
        all_loc = all_Trp;
    case 'Tyr'
        all_loc = all_Tyr;
    case 'Val'
        all_loc = all_Val;
end
if strcmp(which_database,'Dun')
    all_loc = all_loc(:,[1,3]);
else
    all_loc = all_loc(:,1:2);
end

%load('/Users/jennifergaines/Documents/SASA_paper/Dun_1A/Val_locations.mat')
%ind0 = ismember(all_locations(:,2), {'ILE', 'LEU', 'MET', 'PHE', 'VAL', 'THR', 'TRP', 'TYR','SER'});
%uncharged_core = all_locations(ind0,:);
%uncharged_core(:,4:53) = num2cell(1234);
num_samples = 50;
sample_size = 50;
num_babl_sampled = 300;
all_res_dchi = cell(size(all_loc,1), num_samples+2);
all_res_dchi(:,2) = num2cell(0);

for i =  1:size(all_loc,1)
   if mod(i,50)==0
       i
   end
    rng(1);
    
    if exist(strcat(data_folder, lower(all_loc{i,1}), '_', what_res, num2str(all_loc{i,2}),sub_name, '_single_rotation_minE.mat'))
        load(strcat(data_folder, lower(all_loc{i,1}), '_',what_res, num2str(all_loc{i,2}),sub_name, '_single_rotation_minE.mat'));
        if size(unique(all_data(:,1)),1) == num_babl_sampled
            if mod(i,50)==0
                all_data(1,:)
            end
            x = all_data(:,2);
            all_data(:,2:size(all_data,2)-1) = all_data(:,3:size(all_data,2));
            all_data(:,size(all_data,2)) = x;
            all_data(:,size(all_data,2)+1) = exp(-all_data(:,size(all_data,2))*100);
            if ismember(Resi_name , {'ILE' ; 'LEU'; 'PHE'; 'TYR'; 'TRP'})
                
                DOF =2;
            elseif Resi_name == 'MET'
                DOF = 3;
            else
                DOF = 1;
            end
            
            best_dih = [];
            all_dchi = [];
            for randloop = 1:num_samples
                to_sample = randperm(num_babl_sampled,sample_size);
                all_dih = [];
                
                for j = 1:sample_size
                    ind0 = ismember(all_data(:,1), to_sample(j));
                    this_dih = all_data(ind0,2:size(all_data,2));
                    if size(this_dih,1)>0
                        P_E = this_dih(1,size(this_dih,2));
                        if mod(i,50)==0
                            P_E
                        end
                        if j == 1
                            all_dih = this_dih(:,1:DOF);
                            all_dih(:,DOF+1) = P_E;
                        else
                            if size(unique([this_dih(:,1:DOF);all_dih(:,1:DOF)], 'rows'),1) == (size(this_dih,1) + size(all_dih,1))
                                all_dih = [all_dih;this_dih(:,1:DOF), repmat(P_E,size(this_dih,1),1)];
                                
                            else%if size(all_dih,1) >= size(unique_dihedrals,1)
                                
                                
                                ind0 = ismember(all_dih(:,1:DOF), this_dih(:,1:DOF),'rows');
                                all_dih(ind0,DOF+1) = all_dih(ind0,DOF+1)+P_E;
                                ind1 = ismember(this_dih(:,1:DOF), all_dih(:,1:DOF), 'rows');
                                all_dih = [all_dih; this_dih(~ind1,1:DOF), repmat(P_E, sum(~ind1),1)];
                            end
                        end
                    end
                end
                max_val = max(all_dih(:,DOF+1));
                ind0 = ismember(all_dih(:,DOF+1), max_val);
                best_dih = [best_dih;all_dih(ind0,1:DOF),repmat(randloop,sum(ind0),1),  find(ind0 == 1)];
                
                chi1 = all_dih(ind0,1:DOF);
                rand_numb1 = randi(size(chi1,1), 1);
                chi1 = chi1(rand_numb1,:);
                load(strcat(data_folder ,lower(all_loc{i,1}), '_',  what_res, num2str( all_loc{i,2}), '_original.mat'));
                if DOF == 2
                    if (Resi_name == 'PHE' | Resi_name== 'TYR') & orig(2) > 180
                        orig(2) = orig(2)-180;
                    end
                    d1 = abs(chi1(:,1) - orig(1));
                    d2 = abs(chi1(:,2) - orig(2));
                    ind0 = d1 > 180;
                    d1(ind0) = 360-d1(ind0);
                    ind0 = d2 > 180;
                    d2(ind0)  = 360-d2(ind0);
                    this_res = min(sqrt(d1.^2 + d2.^2));
                    
                elseif DOF == 1
                    
                    d1 = abs(chi1(1)-orig(1));
                    ind0 = d1 > 180;
                    d1(ind0) = 360-d1(ind0);
                    this_res = min(d1);
                    
                    
                elseif DOF == 3
                    d1 = abs(chi1(:,1) - orig(1));
                    d2 = abs(chi1(:,2) - orig(2));
                    d3 = abs(chi1(:,3) -orig(3));
                    ind0 = d1 > 180;
                    d1(ind0) = 360-d1(ind0);
                    ind0 = d2 > 180;
                    d2(ind0)  = 360-d2(ind0);
                    ind0 = d3 > 180;
                    d3(ind0)  = 360-d3(ind0);
                    this_res = min(sqrt(d1.^2 + d2.^2 + d3.^2));
                    
                    
                end
                all_dchi = [all_dchi;this_res];
            end
        end
        all_res_dchi(i,1:2) = all_loc(i,1:2);
        all_res_dchi(i,3:num_samples+2) = num2cell(all_dchi');
    end
end
save(strcat(data_folder,what_res, sub_name,'_all_dchi_values.mat'), 'all_res_dchi');
end