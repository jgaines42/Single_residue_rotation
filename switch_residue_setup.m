switch(res_name)
    case 'Ala'
        
        DOF = 1;
        iChi1Array=[4,5,8,11];
        moveAtomID2 = [11:13];
        numAtom = 16;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;2;1;1;1;1;3;5;2];
    case 'Ile'
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,11];
        moveAtomID=[11,15:16,20:22]; %%%%%%%% rotate chi2
        moveAtomID2=[9:11,14:22]; %%%%%%%% rotate chi1
        HG_Array_1=[5,8,10, 17]; %iChi1Array
        HG_Array_2 = [8,9,11, 20];
        moveAtomID_HG1=[17:19];
        moveAtomID_HG2=[20:22];
        numAtom =25;
        DOF =2;
        OH = 0;
        CH3 =2;
        size_index = [5;6;4;3;5;6;4;5;5;5;5;2;1;1;1;1;1;1;1;1;1;1;3;5;2];
        
    case 'Leu'
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,10];
        moveAtomID=[10,11,16:22]; %%%%%%%% rotate chi2
        moveAtomID2=[9:11,14:22]; %%%%%%%% rotate chi1
        HG_Array_1=[8,9,10, 17]; %iHG1Array
        HG_Array_2 = [8,9,11, 20];
        moveAtomID_HG1=[17:19];
        moveAtomID_HG2=[20:22];
        numAtom =25;
        DOF =2;
        OH = 0;
        CH3 =2;
        size_index = [5;6;4;3;5;6;4;5;5;5;5;2;1;1;1;1;1;1;1;1;1;1;3;5;2];
        
    case 'Val'
        iChi1Array = [4,5,8,9];
        moveAtomID2 = [9,10,13,14,15,16,17,18,19];
        HG_Array_1 = [5,8,9,14];
        HG_Array_2 = [5,8,10,17];
        moveAtomID_HG1 = [14,15,16];
        moveAtomID_HG2 = [17:19];
        numAtom = 22;
        DOF = 1;
        OH = 0;
        CH3 = 2;
        size_index = [5;6;4;3;5;6;4;5;5;5;2;1;1;1;1;1;1;1;1;3;5;2];
        
    case 'Phe'
        numAtom = 26;
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,10];
        moveAtomID=[10:14,19:23]; %%%%%%%% rotate chi2
        moveAtomID2=[9:14,17:23]; %%%%%%%% rotate chi1
        DOF = 2;
        OH =0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;10;10;10;10;10;10;2;1;1;1;11;11;11;11;11;3;5;2];
        
    case 'Trp'
        numAtom = 30;
        iChi1Array = [4 5 8 9];
        iChi2Array = [5 8 9 10];
        moveAtomID2 = [9:17, 20:27];
        moveAtomID = [10:17, 22:27];
        DOF = 2;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;10;10;10;16;10;10;10;10;10;2;1;1;1;11;2;11;11;11;11;3;5;2];
        
    case 'Tyr'
        numAtom = 27;
        iChi1Array = [4 5 8 9];
        iChi2Array = [5 8 9 10];
        iOHArray = [12 14 15 24];
        moveAtomID2 = [9:15, 18:24];
        moveAtomID = [10:15, 20:24];
        moveAtomOH = [24];
        DOF = 2;
        OH = 1;
        CH3 = 0;
        size_index = [5; 6; 4 ;3; 5; 6 ;4 ;5; 10;10;10;10;10;10; 12; 2; 1; 1 ;1 ; 11;11;11;11; 2 ;3; 5 ;2];
    case 'Asn'
        iChi1Array=[4,5,8,9];
        iChi2Array=[5,8,9,10];
        moveAtomID=[10,11,16:17]; %%%%%%%% rotate chi2
        moveAtomID2=[9:11,14:17]; %%%%%%%% rotate chi1
        % Ignore terminal H because not a clear distinction of which branch
        % is N and which is O
        
        numAtom =20;
        DOF =2;
        OH = 0;
        CH3 =0;
        size_index = [5; 6; 4 ;3; 5; 6 ;4; 5; 13; 14; 17; 2;1;1;1;2;2;3;5;2];
        
    case 'Cys'
        numAtom = 17;
        iChi1Array = [4, 5, 8, 9];
        iOHArray = [5,8,9,14];
        moveAtomID2 = [9,12,13,14];
        moveAtomOH = 14;
        DOF = 1;
        OH = 1;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;8;2;1;1;1;9;3;5;2];
        
    case 'Glu'
        numAtom = 21;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8,9,10,11];
        moveAtomID = [10:12, 17:18];
        moveAtomID2 = [9:12, 15:18];
        moveAtomID3 = [11,12];
        DOF = 3;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;5;13;14;14;2;1;1;1;1;1;3;5;2];
        
    case 'Met'
        numAtom = 23;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8,9,10,11];
        HG_Array_1 = [9,10,11,18];
        moveAtomID = [10,11,16,17,18,19,20];
        moveAtomID2 = [9,10,11,14,15,16,17,18,19,20];
        moveAtomID3 = [11,18,19,20];
        moveAtomID_HG1 = [18:20];
        DOF = 3;
        OH = 0;
        CH3 = 1;
        size_index = [5;6;4;3;5;6;4;5;5;8;5;2;1;1;1;1;1;1;1;1;3;5;2];
        
    case 'Ser'
        numAtom = 17;
        iChi1Array = [4,5,8,9];
        iOHArray = [5,8,9,14];
        moveAtomID2 = [9,12,13,14];
        moveAtomOH = [14];
        DOF = 1;
        OH = 1;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;15;7;2;1;1;1;16;3;5;2];
        
    case 'Thr'
        numAtom = 20;
        iChi1Array = [4,5,8,9];
        iOHArray = [5,8,9,14];
        HG_Array_1 = [5,8,10,15];
        moveAtomID2 = [9,10,13,14,15,16,17];
        moveAtomOH = [14];
        moveAtomID_HG1 = [15:17];
        DOF = 1;
        OH = 1;
        CH3 = 1;
        size_index = [5;6;4;3;5;6;4;15;7;5;2;1;1;16;1;1;1;3;5;2];
        
    case 'Asp'
        numAtom = 18;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        % iOHArray = [5,8,9,10];
        moveAtomID2 = [9,10,11,14,15];
        moveAtomID = [10,11];
        DOF = 2;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;13;14;14;2;1;1;1;3;5;2];
        
    case 'Gln'
        numAtom = 23;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8,9,10,11];
        moveAtomID = [10:12, 17:10];
        moveAtomID2 = [9:12, 15:20];
        moveAtomID3 = [11,12, 19:20];
        DOF = 3;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;5;13;14;17;2;1;1;1;1;1;2;2;3;5;2];
        
    case 'Arg'
        numAtom = 30;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8,9,10,11];
        iChi4Array = [9,10,11,12];
        iChi5Array = [10,11,12,13];
        moveAtomID = [10:14, 19:27];
        moveAtomID2 = [9:14, 17:27];
        moveAtomID3 = [11:14, 21:27];
        moveAtomID4 = [12:14, 23:27];
        moveAtomID5 = [13:14, 24:27];
        DOF = 5;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;5;5;3;18;17;17;2;1;1;1;1;1;1;1;2;2;2;2;2;3;5;2];
        
    case 'His'
        numAtom = 22;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        moveAtomID2 = [9:13, 16:19];
        moveAtomID = [10:13, 18:19];
        DOF = 2;
        OH = 0;
        CH3 = 0;
        size_index = [5;6;4;3;5;6;4;5;10;16;10;10;16;2;1;1;1;11;11;3;5;2];
        
    case 'Lys'
        numAtom = 28;
        iChi1Array = [4,5,8,9];
        iChi2Array = [5,8,9,10];
        iChi3Array = [8 9 10 11];
        iChi4Array = [9 10 11 12];
        HG_Array_1 = [10 11 12 23];
        moveAtomID2 = [9:12, 15:25];
        moveAtomID = [10:12, 17:25];
        moveAtomID3 = [11:12, 19:25];
        moveAtomID4 = [12, 21:25];
        moveAtomID_HG1 = [23:25];
        DOF = 4;
        OH = 0;
        CH3 = 1;
        size_index = [5;6;4;3;5;6;4;5;5;5;5;3;2;1;1;1;1;1;1;1;1;1;2;2;2;3;5;2];
        
    case 'Gly'
        fprintf('Not yet supported\n' );
    case 'Pro'
        fprintf('Not yet supported\n' );
    otherwise
        fprintf('Invalid amino acid\n' );
end
