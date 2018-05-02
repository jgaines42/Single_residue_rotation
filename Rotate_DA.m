%Rotate_DA
function new_Pos = Rotate_DA(Pos, ang, SA, dt, iCA, mAID)

    Position=Pos; % Coordinates of the dipeptide
    setChi = ang; % Angle to rotate them by
    delta_term = dt; % dt= pi*sign(InitChi)*InitChi1/180
    iChiArray = iCA; % indexes of the atoms that define the dihedral angle
    moveAtomID = mAID; % Atoms to rotate 
    subtract_array = SA; % SA = repmat(Position(iChiArray(2),:),numAtom,1)
    
    deltaChi1_F153=delta_term-setChi*pi/180;

    TempPosition=Position-subtract_array;
    CAtoCB_F153=-TempPosition(iChiArray(3),:);
    CAtoCB_F153=CAtoCB_F153./norm(CAtoCB_F153);


    q0 = cos(deltaChi1_F153/2);
    q1 = CAtoCB_F153(1)*sin(deltaChi1_F153/2);
    q2 = CAtoCB_F153(2)*sin(deltaChi1_F153/2);
    q3 = CAtoCB_F153(3)*sin(deltaChi1_F153/2);
    Q =    [(q0^2 + q1^2 - q2^2 - q3^2), 2*(q1*q2 - q0*q3),           2*(q0*q2 + q1*q3);
        2*(q1*q2 + q0*q3),           (q0^2 - q1^2 + q2^2 - q3^2), 2*(-q0*q1 + q2*q3);
        2*(-q0*q2 + q1*q3),          2*(q0*q1 + q2*q3), (q0^2 - q1^2 - q2^2 + q3^2)];
    newPos=Q*TempPosition(moveAtomID(:),:)';
    TempPosition(moveAtomID(:),:)=newPos';
    new_Pos=TempPosition+subtract_array;
end