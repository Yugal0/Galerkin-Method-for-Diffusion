clc
clear

%The Experimental Data (D[m^2/s] vs Li Saturation in NMC811 Cathode
%Material) from Ball et. al.
C_Exp=0:0.01:1;
D_Exp=[1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.97000000000000e-15	1.95871794871795e-15	1.93246153846154e-15	1.88057988240542e-15	1.82604733885277e-15	1.77151479530013e-15	1.75000000000000e-15	1.75000000000000e-15	1.75000000000000e-15	1.75000000000000e-15	1.75000000000000e-15	1.75000000000000e-15	1.72415384688095e-15	1.69039560607239e-15	1.67443076940527e-15	1.66000000000000e-15	1.66000000000000e-15	1.65677948731521e-15	1.64627692279901e-15	1.64000000000000e-15	1.64000000000000e-15	1.62076923231769e-15	1.54199999982496e-15	1.46323076733223e-15	1.36872435659487e-15	1.27354486941539e-15	1.17836538223590e-15	1.09518090133021e-15	1.01315225578639e-15	9.31123610242575e-16	8.60628439783704e-16	8.28643358007756e-16	7.96658276231807e-16	7.64673194455859e-16	7.33874938039106e-16	7.30825806467350e-16	7.27776674895595e-16	7.24727543323840e-16	7.15605591488646e-16	6.54977154345584e-16	5.94348717202522e-16	5.33720280059460e-16	4.73091842916398e-16	4.45872820826228e-16	4.21716923073207e-16	3.97561025320186e-16	3.70044377402405e-16	3.14602958123884e-16	2.59161538845363e-16	2.03720119566843e-16	1.66534850545147e-16	1.33725753836087e-16	1.00916657127028e-16	8.67230773551425e-17	7.49076925607793e-17	6.30923077664160e-17	5.12769229720528e-17	3.94615381776896e-17	2.76461533833264e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17	1.92000000000000e-17];



%For Fitting Poly Diffusisivty in Exp data was noramlized by 2e-14
function [Fval]=givePolyFit1(C)
p1  =  -15.8080;  
p2  =  61.2357;  
p3   = -97.2872;  
p4   = 82.0729 ;   
p5   = -38.9613    ;
p6   = 9.7527   ;    
p7   = -1.1536  ;  
p8   = 0.0518  ;    
p9   = 0.0980 ;      

x=C;
%Fitted Poly
Fval=p1*x.^8 + p2*x.^7 + p3*x.^6 + p4*x.^5 + p5*x.^4 + p6*x.^3 + p7*x.^2 + p8.*x + p9;
end


function [Fval]=givePolyFit2(C)
p1   = 0.1800 ;   
p2  =  -3.7603   ; 
p3  =  10.1035 ;   
p4  =  -9.9922  ;  
p5  =  4.0525  ;  
p6  =  -0.7328 ;   
p7  =  0.0521  ;   
p8  =  0.0977 ;   

x=C;
%Fitted Poly
Fval=p1*x.^7 + p2*x.^6 + p3*x.^5 + p4*x.^4 + p5*x.^3 + p6*x.^2 + p7*x + p8;

end


function [Fval]=givePolyFit3(C)

p1=-3.1302;
p2=9.2325;
p3=-9.3898;
p4=3.8353;
p5=-0.6944;
p6=0.0494;
p7=0.0977;

%Fitted Poly
Fval=p1*C.^6 + p2*C.^5 + p3*C.^4 + p4*C.^3 + p5*C.^2 + p6*C + p7;
end


function [Fval]=givePolyFit4(C)
p1=0.7936;
p2=-1.3360;
p3=0.5023;
p4=-0.0569;
p5=0.0998;


%Fitted Poly
Fval=p1*C.^4 + p2*C.^3 + p3*C.^2 + p4*C.^1 + p5;
end


function [Fval]=givePolyFit5(C)

p1=-0.1289;
p2=-0.13;

%Fitted Poly
Fval=p1*C.^4 - p2;
end


function [Fval]=giveBestFit(C)
a=0.1008;
b=-8.2105;
c=0.6986 ;
d=0.0009;
Fval = d + (a-d)*exp(-exp(-b*(C-c)));
end

% Taylor Series Poly from Best Fit Exponential Expression

% syms CC 
% a=0.1008;
% b=-8.2105;
% c=0.6986 ;
% d=0.0009;
% Fval = d + (a-d)*exp(-exp(-b*(CC-c)));
% givePolyFitfromBest=matlabFunction(vpa(taylor(Fval,CC,0,"Order",10)))



%Plot Fitted Poly vs Experiment 
figure 
subplot(1,2,1)

%The Conc. to plot at
C=C_Exp;

plot(C,givePolyFit1(C),"b")
hold on
plot(C,givePolyFit2(C),"-.m")
hold on
plot(C,givePolyFit3(C),"-.c")
hold on
plot(C,givePolyFit4(C),"-.g")
hold on
plot(C,givePolyFit5(C),":b")
hold on
plot(C,giveBestFit(C),"--r")
hold on
% plot(C,givePolyFitfromBest(C),"r")
% hold on
plot(C,D_Exp/2e-14,".k")
xlabel("C")
ylabel("f")
% legend(["Fitted Poly 8 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Best Fit (But has exponential)", "Poly Fit from Taylor Series of Best Fit","Experiment"])
% legend(["Fitted Poly 8 Order","Fitted Poly 7 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Best Fit (But has exponential)", "Experiment"])
legend(["Fitted Poly 8 Order","Fitted Poly 7 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Simple Poly 4th Order","Best Fit (But has exponential)", "Experiment"])
title("Normal Scale")

subplot(1,2,2)
semilogy(C,givePolyFit1(C),"b")
hold on
semilogy(C,givePolyFit2(C),"-.m")
hold on
semilogy(C,givePolyFit3(C),"-.c")
hold on
semilogy(C,givePolyFit4(C),"-.g")
hold on
semilogy(C,givePolyFit5(C),":b")
hold on
semilogy(C,giveBestFit(C),"--r")
hold on
% semilogy(C,givePolyFitfromBest(C),"r")
% hold on
semilogy(C,D_Exp/2e-14,".k")
hold on
yline(1e-1,"--")
hold on
yline(1e-3,"--")
xlabel("C")
ylabel("f")
% legend(["Fitted Poly 8 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Best Fit (But has exponential)", "Poly Fit from Taylor Series of Best Fit","Experiment"])
% legend(["Fitted Poly 8 Order","Fitted Poly 7 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Best Fit (But has exponential)", "Experiment"])
legend(["Fitted Poly 8 Order","Fitted Poly 7 Order","Fitted Poly 6 Order","Fitted Poly 4 Order","Simple Poly 4th Order","Best Fit (But has exponential)", "Experiment"])
title("Log Scale in Y-axis")

%% 

figure;
semilogy(C,D_Exp,".k")
xlabel("C")
ylabel("f")
