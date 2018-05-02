
if DOF >= 3
    load(strcat('chi3_', res_name, '_clash.mat'));
    c3_clashList = n;
    if next_pro == 1
        ind0 = c3_clashList(:,2) == numAtom;
        c3_clashList = c3_clashList(~ind0,:);
    end
    
    c3_Clash = [c3_clashList, (Atom_sizes(c3_clashList(:,1),:)+Atom_sizes(c3_clashList(:,2),:))];
    InitChi3=calcDA2(iChi3Array,Position); %%%%%%%%
    InitChi3=mod(real(InitChi3),360);
    
    
    
    ind0 = ismember(moveAtomID3, n);
    moves_here = moveAtomID3(ind0);
    c3_Clash(:,3)=c3_Clash(:,3).^2;
end

if DOF >= 2
    load(strcat('chi2_', res_name, '_clash.mat'));
    c2_clashList = n;
    if next_pro == 1
        ind0 = (c2_clashList(:,2) == numAtom);
        c2_clashList = c2_clashList(~ind0,:);
    end
    c2_Clash = [c2_clashList, (Atom_sizes(c2_clashList(:,1),:)+Atom_sizes(c2_clashList(:,2),:))];
    c2_Clash(:,3)=c2_Clash(:,3).^2;
    
    
    InitChi2=calcDA2(iChi2Array,Position); %%%%%%%%
    InitChi2=mod(real(InitChi2),360);
    max_Chi2 = 72;
    if strcmp(res_name,'Phe') || strcmp(res_name, 'Tyr') || strcmp(res_name, 'Asp')
        max_Chi2 = 72/2;
    end
    
    ind0 = ismember(moveAtomID, n);
    moves_here = moveAtomID(ind0);
    
    
end
if DOF >= 1
    load(strcat('chi1_', res_name, '_clash.mat'));
    c1_clashList = n;
    if next_pro == 1
        ind0 = c1_clashList(:,2) == numAtom;
        c1_clashList = c1_clashList(~ind0,:);
    end
    c1_Clash = [c1_clashList, (Atom_sizes(c1_clashList(:,1),:)+Atom_sizes(c1_clashList(:,2),:))];
    
    
    InitChi1=calcDA2(iChi1Array,Position);
    InitChi1=mod(real(InitChi1),360);
    
    ind0 = ismember(moveAtomID2, n);
    moves_here = moveAtomID2(ind0);
    c1_Clash(:,3) = c1_Clash(:,3).^2;
    
end
if CH3 == 2
    load(strcat('HD2_', res_name, '_clash.mat'));
    HG2_clashList = n;
    if next_pro == 1
        ind0 = HG2_clashList(:,2) == numAtom;
        HG2_clashList = HG2_clashList(~ind0,:);
    end
    HG2_Clash = [HG2_clashList, (Atom_sizes(HG2_clashList(:,1),:)+Atom_sizes(HG2_clashList(:,2),:))];
    InitChi_HG2=calcDA2(HG_Array_2,Position);
    InitChi_HG2=mod(real(InitChi_HG2),360);
    
    
    
    ind0 = ismember(moveAtomID_HG2, n);
    moves_here = moveAtomID_HG2(ind0);
    HG2_Clash(:,3) = HG2_Clash(:,3).^2;
    
end
if CH3 >=1
    load(strcat('HD1_', res_name, '_clash.mat'));
    HG1_clashList = n;
    if next_pro == 1
        ind0 = HG1_clashList(:,2) == numAtom;
        HG1_clashList = HG1_clashList(~ind0,:);
    end
    HG1_Clash = [HG1_clashList, (Atom_sizes(HG1_clashList(:,1),:)+Atom_sizes(HG1_clashList(:,2),:))];
    InitChi_HG1=calcDA2(HG_Array_1,Position); %InitChi1
    InitChi_HG1=mod(real(InitChi_HG1),360);
    
    
    
    
    ind0 = ismember(moveAtomID_HG1, n);
    moves_here = moveAtomID_HG1(ind0);
    HG1_Clash(:,3) = HG1_Clash(:,3).^2;
end
if OH == 1
    if res_name == 'Cys'
        load(strcat('SH_', res_name, '_clash.mat'));
    else
        load(strcat('OH_', res_name, '_clash.mat'));
    end
    OH_clashList = n;
    if next_pro == 1
        ind0 = OH_clashList(:,2) == numAtom;
        OH_clashList = OH_clashList(~ind0,:);
    end
    OH_Clash = [OH_clashList, (Atom_sizes(OH_clashList(:,1),:)+Atom_sizes(OH_clashList(:,2),:))];
    InitOH=calcDA2(iOHArray,Position); %InitChi1
    InitOH=mod(real(InitOH),360);
    
    
    
    
    ind0 = ismember(moveAtomOH, n);
    moves_here = moveAtomOH(ind0);
    OH_Clash(:,3) = OH_Clash(:,3).^2;
end
