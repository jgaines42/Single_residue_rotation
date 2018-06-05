%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function new_Pos = Rotate_DA(Pos, ang, SA, dt, iCA, mAID)
% 
% Rotates the coordinates by the dihedral angle given
%
% Input:
% Position: Coordinates of the dipeptide
% SetChi:  Angle to rotate to
% subtract_array: array used to move second atom of dihedral angle to the
% orgin
% delta_term: pi*sign(InitChi)*InitChi/180
% iChiArray: indexes of the atoms that define the dihedral angle
% moveAtomID: Atoms to move during rotation
%
% Return:
% new_Pos: new coordinates of dipeptide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function new_Pos = Rotate_DA(Position, setChi, subtract_array, delta_term, iChiArray, moveAtomID)
    
    deltaChi1_F153=delta_term-setChi*pi/180;

    TempPosition=Position-subtract_array; % Move to orgin
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
    new_Pos=TempPosition+subtract_array; %Move back to original location
end