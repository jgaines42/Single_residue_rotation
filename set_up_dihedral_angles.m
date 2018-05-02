%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sampledCombin = set_up_dihedral_angles(DOF)
%
% Sets up the storage array for each residue
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sampledCombin = set_up_dihedral_angles(DOF)
        if DOF == 3 || DOF == 4
            chi1Vary=[5:5:360];
            list1=chi1Vary';
            sampledCombin= zeros(72*72*72, 4);
            for i=1: size(list1,1)
                for j=1: size(list1,1)
                    for k=1: size(list1,1)
                        sampledCombin(k +72*(j-1) + 72*72*(i-1), 1) = list1(i);
                        sampledCombin(k +72*(j-1) + 72*72*(i-1), 2) = list1(j);
                        sampledCombin(k +72*(j-1) + 72*72*(i-1), 3) = list1(k);
                    end
                end
            end
            if DOF == 4
                sampledCombin = [repmat(sampledCombin,72,1), zeros(72^4,1)];
                for i = 1:size(list1,1)
                    sampledCombin((i-1)*72 + 1 : i*72,4) = i;
                end
            end
            
        elseif DOF == 2
            chi1Vary=5:5:360;
            sampledCombin= zeros(72*72, 4);
            for i = 1:72
                sampledCombin((i-1)*72+1:i*72,1) = chi1Vary(i);
                sampledCombin((i-1)*72+1:i*72,2) = chi1Vary(:);
            end
        else
            
            chi1Vary=5:5:360;
            sampledCombin= [chi1Vary', zeros(72,2)];
        end
end