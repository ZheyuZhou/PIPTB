function [V,T] = autoGen_PIPTBEnergy(in1,in2,uw,mK,IK,rK,mA,IA,lAC,mL,g)
%AUTOGEN_PIPTBENERGY
%    [V,T] = AUTOGEN_PIPTBENERGY(IN1,IN2,UW,MK,IK,RK,MA,IA,LAC,ML,G)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Apr-2021 13:22:42

dphi = in2(2,:);
dtheta = in2(1,:);
theta = in1(1,:);
t2 = cos(theta);
V = g.*lAC.*mA.*t2;
if nargout > 1
    t3 = dphi.^2;
    t4 = dtheta.^2;
    t5 = rK.^2;
    T = (IA.*t4)./2.0+(IK.*t3)./2.0+(mA.*((dphi.*rK+dtheta.*lAC.*t2).^2+lAC.^2.*t4.*sin(theta).^2))./2.0+(mK.*t3.*t5)./2.0+(mL.*t3.*t5)./2.0;
end
