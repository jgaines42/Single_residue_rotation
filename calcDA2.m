function [DA]=calcDA2(AtomArray, position)
	% disp([i1, i2, i3, i4])
	rijxrkj2=0.0;rjkxrlk2=0.0;volijkl=0.0;
	rijxrkjxrjkxrlk2=0.0;
	
	i1=AtomArray(1);
	i2=AtomArray(2);
	i3=AtomArray(3);
	i4=AtomArray(4);
	 
	subij=position(i1,:)-position(i2,:);
	subjk=position(i2,:)-position(i3,:);
	subkl=position(i3,:)-position(i4,:);
	%subij = b1;
	%subjk = b2;
	%subkl = b3;
	
	rijxrkj= cross(subij,subjk);
	rjkxrlk= cross(subjk,subkl);
	rijxrkj2=sum(rijxrkj.^2); %|rij x rkj|^2
	rjkxrlk2=sum(rjkxrlk.^2); %|rjk x rlk|^2
	volijkl=  	rijxrkj(1)*rjkxrlk(1)+...
				rijxrkj(2)*rjkxrlk(2)+...
				rijxrkj(3)*rjkxrlk(3);%volijkl=(rij x rkj)*(rjk x rlk),for calc cosphi
		
	rkj=sqrt(sum(subjk.^2));
	
	rijxrkjxrjkxrlk=cross(rijxrkj,rjkxrlk);
	cosphi=volijkl/sqrt(rijxrkj2*rjkxrlk2);
	sinphi=  	(rijxrkjxrjkxrlk(1)*subjk(1)+...
				rijxrkjxrjkxrlk(2)*subjk(2)+...
				rijxrkjxrjkxrlk(3)*subjk(3))/sqrt(rijxrkj2*rjkxrlk2)/rkj;
	
	phi=-sign(sinphi)*acos(cosphi);
	
	DA=phi/pi*180; 
end
