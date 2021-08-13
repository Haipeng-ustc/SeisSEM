%
% Compute the components of a moment tensor for a fault mechanism
%
% INPUT
%
% str	:	strike of a fault mechanism
% dip	:	dip of a fault mechanism
% rake	:	rake of a fault mechanism
%
% OUTPUT
%
% CM_FD :	the moment tensor component for finite difference modeling
%               You should directly paste the CM_FD in the corresponding part
%		of fd3dinput.asc to define the moment tensor source
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CM] = Moment_Tensor_Solution(strike, dip, rake)

pi180 = pi/180;
CS=cos(strike*pi180);  %calculate each moment component
SS=sin(strike*pi180);
CDI=cos(dip*pi180);
SDI=sin(dip*pi180);
CR=cos(rake*pi180);
SR=sin(rake*pi180);
AS1=CR*CS+SR*CDI*SS;
AS2=CR*SS-SR*CDI*CS;
AS3=-SR*SDI;
AN1=-SDI*SS;
AN2=SDI*CS;
AN3=-CDI;
CM11=2.*AS1*AN1;  %change to the normal convention by eliminating the minus sign
CM22=2.*AS2*AN2;
CM33=2.*AS3*AN3;
CM12=(AS1*AN2+AS2*AN1);
CM13=(AS1*AN3+AS3*AN1);
CM23=(AS2*AN3+AS3*AN2);

CM = [ CM11 CM12 CM13;
          CM12 CM22 CM23;
          CM13 CM23 CM33];

end 