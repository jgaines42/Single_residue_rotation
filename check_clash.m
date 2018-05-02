function [energy] = check_clash(prior_residues, new_residue)
    size_prior = size(prior_residues,1);
    size_new = size(new_residue,1);
    all_ind = repmat([1:size_prior]',size_new,1);
    nn = repmat([1:size_new], size_prior,1);
    all_ind(:,2) = reshape(nn,[size(all_ind,1),1]);
  
    dist = prior_residues(all_ind(:,1),1:3)-new_residue(all_ind(:,2),1:3);
    distemp = sum(dist.^2,2);
    DA_Clash = (prior_residues(all_ind(:,1),4)+new_residue(all_ind(:,2), 4)).^2;
    clash = distemp - DA_Clash;
    ind_clash = find(clash < 0);
    if size(ind_clash,1) > 0
        energy = sum((1 - (DA_Clash(ind_clash,1)./distemp(ind_clash)).^3).^2);
    else
        energy = 0;
    end
end