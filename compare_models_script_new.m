%% Creating Variables

%clear memory and command window
clear all
clc
dd=8%digit precision
digits(dd)

% Approximating with first N ortho polynomials
N=3;  %minimum is N=2, max allowed is N=6
tic %timer on part 1

% Create symbolic vars for r and t
syms r t
assume(r,"real")
assume(t,"real")

% P will be the array of ortho polys
p=sym("P",[1 N]);
P= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p, 'UniformOutput', false));


% alpha will be the array of coeff of ortho poly
aa=sym("alpha",[1 N]);
alpha= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa, 'UniformOutput', false));


% Approximation of conc. inside particle C(r,t) with ortho polys
Capprox = sym("0");  
for ii=1:N    
    Capprox=Capprox+P(ii)*alpha(ii);
end

% See the Capprox
Capprox;

%initializing with zero
polys_to_use=sym("0");

%array to contain the unknown coeff of ortho polys
poly_coeff=sym("nu",[N N]);

%creating template of ortho polys with still unknown coeff
for ii=1:N
    if(ii==1)
        polys_to_use(ii)=1;
    else
        polys_to_use(ii)=sym("0");
        for jj=2:ii
            polys_to_use(ii)=polys_to_use(ii)+poly_coeff(ii,jj)*r^jj;
        end
        polys_to_use(ii)=1+polys_to_use(ii);
    end
end
%see ortho poly template
polys_to_use';

kk=1;
poly_coeff_in_use=sym("0");
for ii=2:N
    for jj=2:ii
        poly_coeff_in_use(kk)=poly_coeff(ii,jj);
        kk=kk+1;
    end
end

% Orthogonality Conditions (which will help in determining the unknown coeff of ortho polys)
ortho_eqn=sym("0");
kk=1;
for ii=1:N
    for jj=1:N
        if(ii~=jj && ii<jj)           %lower half of full matrix
            ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
            kk=kk+1;
        end
    end
end


%Need to assume coeff to be real here, coz somehow MATLAB forgets if I
%assume them real just on definition
assume(poly_coeff,"real");

%solve the orthogonality equations to get the unknown coeff of ortho polys
poly_coeff_sol=solve(ortho_eqn);

%something to make output of solve command more usable (the output is a MATLAB structure by default, but we want to work with MATLAB sym array (symbolic array))
if(N>2)  % when N=2 we have just one unkown and hence somehow the solve command returns a number instead of a structure (containing multiple numbers)
    poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';

end
% poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';

%Substitute the values of unknown coeff in the ortho polys expressions
polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);

JP_polys_to_use=polys_to_use;
%see the final ortho polys to be used 
JP_polys_to_use';

%Write Capprox in explicit form by replacing P1(r) P2(r)... polys with
%their expression in terms of r
Capprox_JP=subs(Capprox,P,JP_polys_to_use);
% Generate Symmetric Like Polynomials

%initializing with zero
polys_to_use=sym("0");

%array to contain the unknown coeff of ortho polys
poly_coeff=sym("nu",[N N]);

%creating template of ortho polys with still unknown coeff
for ii=1:N
    if(ii==1)
        polys_to_use(ii)=1;
    else
        polys_to_use(ii)=sym("0");
        for jj=1:ii-1
            polys_to_use(ii)=polys_to_use(ii)+poly_coeff(ii,jj)*r^(2*jj);
        end
        polys_to_use(ii)=1+polys_to_use(ii);
    end
end
%see ortho poly template
polys_to_use';

%Out of all the coeff in unknown coeff 2D array, only few will have
%non-zero terms (following the pattern as per code below ), Note only those
%with non-zero coeff
kk=1;
poly_coeff_in_use=sym("0");
for ii=2:N
    for jj=1:ii-1
        poly_coeff_in_use(kk)=poly_coeff(ii,jj);
        kk=kk+1;
    end
end

% Orthogonality Conditions (which will help in determining the unknown coeff of ortho polys)
ortho_eqn=sym("0");
kk=1;
for ii=1:N
    for jj=1:N
        if(ii~=jj && ii<jj)           %lower half of full matrix
            ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
            kk=kk+1;
        end
    end
end


%Need to assume coeff to be real here, coz somehow MATLAB forgets if I
%assume them real just on definition
assume(poly_coeff,"real");

%solve the orthogonality equations to get the unknown coeff of ortho polys
poly_coeff_sol=solve(ortho_eqn);

if(N>2)  % when N=2 we have just one unkown and hence somehow the solve command returns a number instead of a structure (containing multiple numbers)
    poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';

end
%something to make output of solve command more usable (the output is a MATLAB structure by default, but we want to work with MATLAB sym array (symbolic array))
%poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';

%Substitute the values of unknown coeff in the ortho polys expressions
polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);

SP_polys_to_use=polys_to_use;

%see the final ortho polys to be used 
SP_polys_to_use';

%Write Capprox in explicit form by replacing P1(r) P2(r)... polys with
%their expression in terms of r
Capprox_SP=subs(Capprox,P,SP_polys_to_use);

% Capprox_SP;
Coeff_SP_poly = coeffs(Capprox_SP,r);
Coeff_SP_poly = subs(Coeff_SP_poly,alpha,aa);
% Generate Thiag Polynomials

Nc=N-1;
thiag_polys_to_use=SP_polys_to_use(2:end);

%see the final ortho polys to be used 
thiag_polys_to_use';
rootP=vpasolve(thiag_polys_to_use(end));
rootP=sort(abs(rootP(1:end/2)));
rootP=[rootP' 1];

%weights in the lagrange polynomial
lgr_wts=sym("phi",[1 Nc+1]);
Lgr_wts= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), lgr_wts, 'UniformOutput', false));

for ii=1:length(rootP)
    Lgr_wts(ii)=1;
    for jj=1:length(rootP)
        if ~(ii==jj)
            Lgr_wts(ii)=Lgr_wts(ii)*(r^2-rootP(jj)^2)/(rootP(ii)^2-rootP(jj)^2);
        end
    end
end
% c will be the array of coeff of ortho poly
cc=sym("c",[1 Nc+1]);
c= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), cc, 'UniformOutput', false));

%Write Capprox in explicit form by replacing P1(r) P2(r)... polys with
%their expression in terms of r
Capprox_thiag=subs(Capprox,[P alpha],[Lgr_wts c]);

time_part_1=toc;

%% Model Parameters and Simulation Parameters
tic
%Parameters for the Model 
%And also the Simulation parameters, Had to write them
%here, But they are not used as much till the final ODE system simulation)
pStruct.FvsCprofile=char("simplePoly");
pStruct.Deltavstauprofile=char("testCrate");
pStruct.C0vsrprofile=char("const");
pStruct.DeltaIsZeroForSeconds=20;
pStruct.Compress=1; 
pStruct.R=2e-6;
pStruct.c0=51554;
pStruct.D0=2e-14;
pStruct.D0multiplier=1e-1;
pStruct.j0=4.3e-6; 
pStruct.CrateMultiplier=3;
pStruct.CrateMultiplierNegative=0.7;
pStruct.C0=1;
pStruct.RelTol=1e-9;
pStruct.AbsTol=1e-9;
pStruct.MaxStep=1e-2;

pStruct.rNodes=1000;
pStruct.tNodes=2000;
pStruct.tauMax=2;
pStruct.totalTimeSeconds=pStruct.tauMax*1000;
rpoints=linspace(0,1,pStruct.rNodes);
pStruct.rpoints=rpoints;
tpoints=linspace(0,pStruct.totalTimeSeconds,pStruct.tNodes);
taupoints=(pStruct.D0/pStruct.R^2)*tpoints;
pStruct.taupoints=taupoints;
%% Applying the Galerkin Method
% Setting up the Variables and Parameters
%% ff is the non-dimensionalized diffusivity number (as per Thiagarajan et. al. paper), ff is a function of conc

ff = @(r,Capprox) firstOutputOnly(r,Capprox,pStruct);

%giveFforC gives two outputs i.e. value Diffusivity and Derivative of
%Diffusivity wrt C But we only want the value of Diffusivity
function out1 = firstOutputOnly(r,Capprox,pStruct)
    [out1, ~] = giveFforC(Capprox, pStruct.FvsCprofile, pStruct);
end
%% Info for the Quadrature

%number of points in quadrature should depend on the highest degree of r (or x) in
%the integrand
% N_quad_JP=ceil((N+1)/2);
% N_quad_SP=4*ceil((2*(N-1)+1)/2);
% % N_quad_SP=50;
% N_quad_thiag=Nc;
% N_quad_thiag_GJ=4*ceil((2*(Nc)+1)/2);
% N_quad_thiag_GJ=6;
% N_quad_thiag_GJ=9;


% New Quad Point Rules

% N_quad_JP=ceil((N+1)/2);
% N_quad_SP=10*N-9;
% N_quad_thiag=Nc;
% N_quad_thiag_GJ=10*N-9;

N_quad_JP=ceil((N+1)/2);
% N_quad_SP=10*N-9;
m=8;
% N_quad_SP=(1+m)*(N-1)+1;
% N_quad_SP=(m+2)*(N-1);
N_quad_SP=8;
N_quad_thiag=Nc;
% N_quad_thiag_GJ=(m+2)*(N-1);
N_quad_thiag_GJ=8;

% N_quad_thiag_GJ=(1+m)*(N-1)+1;


% a and b are the alpha and beta as in standard jacobi polys
a=0;
b=2;

%symbolic function for flux
syms delta(t)

%The weights used in Method of Weighted Residual are set same as the ortho
%polys i.e. As per the Galerkin Method
W_JP=JP_polys_to_use;
W_SP=SP_polys_to_use;
W_thiag=Lgr_wts;
% Galerkin Method to get _N_ number of Model equations 

%Below is the code for the galerkin method equations (the 3 on both sides can be removed now?)
for ii=1:N
    Galerk_eqn_JP(ii)=int(W_JP(ii)*3*r^2*diff(Capprox_JP,t),r,0,1)==subs(W_JP(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_JP(ii),r)*3*ff(r,Capprox_JP)*diff(Capprox_JP,r),N_quad_JP,a,b);
    Galerk_eqn_SP(ii)=int(W_SP(ii)*3*r^2*diff(Capprox_SP,t),r,0,1)==subs(W_SP(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_SP(ii),r)*3*ff(r,Capprox_SP)*diff(Capprox_SP,r),N_quad_SP,a,b);
    Galerk_eqn_thiag(ii)=int(W_thiag(ii)*3*r^2*diff(Capprox_thiag,t),r,0,1)==subs(W_thiag(ii)*3*r^2*(-delta),r,1)-giveThiagQuadResult(diff(W_thiag(ii),r)*ff(r,Capprox_thiag)*diff(Capprox_thiag,r),N_quad_thiag);
    Galerk_eqn_thiag_GJ(ii)=int(W_thiag(ii)*3*r^2*diff(Capprox_thiag,t),r,0,1)==subs(W_thiag(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_thiag(ii),r)*3*ff(r,Capprox_thiag)*diff(Capprox_thiag,r),N_quad_thiag_GJ,a,b);
end

% Slow Manifold 2 Equation Model
SM_2eqn(1)=diff(c(1),t)==-3*delta;
SM_2eqn(2)=diff(c(2),t)==20.5551*c(1)-20.5551*c(2)-7.11102*delta-0.0825422*diff(delta,t)+0.00049392*diff(delta,t,2)-3.01183e-6*diff(delta,t,3);



% Functions written for Quadrature

%function implementing the quadrature, but it needs the points and weights
%to be provided by one other function (which is the next function)
function result = giveQuadResult(integrand,n,a,b)
    syms r
    [quad_pts,quad_wts]=giveGaussJacobiQuad(n,a,b);
    for jj=1:length(quad_pts)
        term_quad(jj)=quad_wts(jj)*subs(integrand,r,quad_pts(jj));
    end
    result=vpa(sum(term_quad));
end

%Provie the points and weights as per gauss Jacobi Quadrature
function [quad_pts,quad_wts]=giveGaussJacobiQuad(n,a,b)
    %need to define symbolic var r again, so that can be used inside
    %function
    syms r

    %jacobiP(n,a,b,r) is MATLAB Builtin command to generate Jacobi poly
    jacobi_poly_n=jacobiP(n,a,b,r);
    % assignin('base',"jacobi_poly_n",jacobi_poly_n);
    %points for the quadrature are the roots of  
    quad_pts_0=vpasolve(jacobi_poly_n);
    % quad_pts_0    
    quad_pts=(quad_pts_0+1)/2; 

    %Weights for Gauss-Jacobi (GJ) Quadrature (from Gauss-Jacobi Quadrature Formula on
    %WiKipedia)
    %the n+1 jacobi poly will be used the weights formula
    jacobi_poly_n_plus_1=jacobiP(n+1,a,b,r);
    % assignin('base',"jacobi_poly_n_plus_1",jacobi_poly_n_plus_1);
    for ii=1:length(quad_pts)
        %the quad points used in the GJ quad formula for quad weights, are the roots for Jacobi polys
        % in (-1,1) (not in (0,1) i.e. not the ones we get after r->(r+1)/2)
        % subs(diff(jacobi_poly_n,r)*jacobi_poly_n_plus_1,r,quad_pts_0(ii))
        quad_wts(ii)=-((2*n+a+b+2)/(n+a+b+1))*((gamma(n+a+1)*gamma(n+b+1))/(gamma(n+a+b+1)*factorial(n+1)))*(2^(a+b))/subs(diff(jacobi_poly_n,r)*jacobi_poly_n_plus_1,r,quad_pts_0(ii));
    end
    

    %Need to compensate for the two substitions we make in original
    %jacobiP(n,0,2,r) polys, the two substitutions are [1] r->(2r-1), which is to make polys
    %orthogonal to r^2 instead of (1+r)^2, [2] r-> (r+1)/2, which is to
    %bring poly from (-1,1) to (0,1)
    % digits(6);
    quad_wts=vpa(quad_wts/(2*4))';
end

function result = giveThiagQuadResult(integrand,Nc)
    syms r
    [quad_pts,quad_wts]=giveThiagQuad(Nc);
    for jj=1:length(quad_pts)
        term_quad(jj)=quad_wts(jj)*subs(integrand,r,quad_pts(jj));
    end
    result=vpa(sum(term_quad));
end


function [quad_pts,quad_wts]=giveThiagQuad(Nc)
    N=Nc+1;

    r=sym("r");

    %initializing with zero
    polys_to_use=sym("0");
    
    %array to contain the unknown coeff of ortho polys
    poly_coeff=sym("nu",[N N]);
    
    %creating template of ortho polys with still unknown coeff
    for ii=1:N
        if(ii==1)
            polys_to_use(ii)=1;
        else
            polys_to_use(ii)=sym("0");
            for jj=1:ii-1
                polys_to_use(ii)=polys_to_use(ii)+poly_coeff(ii,jj)*r^(2*jj);
            end
            polys_to_use(ii)=1+polys_to_use(ii);
        end
    end
    %see ortho poly template
    % polys_to_use'
    %Out of all the coeff in unknown coeff 2D array, only few will have
    %non-zero terms (following the pattern as per code below ), Note only those
    %with non-zero coeff
    kk=1;
    poly_coeff_in_use=sym("0");
    for ii=2:N
        for jj=1:ii-1
            poly_coeff_in_use(kk)=poly_coeff(ii,jj);
            kk=kk+1;
        end
    end
    
    % Orthogonality Conditions (which will help in determining the unknown coeff of ortho polys)
    ortho_eqn=sym("0");
    kk=1;
    for ii=1:N
        for jj=1:N
            if(ii~=jj && ii<jj)           %lower half of full matrix
                ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
                kk=kk+1;
            end
        end
    end
    
    
    %Need to assume coeff to be real here, coz somehow MATLAB forgets if I
    %assume them real just on definition
    assume(poly_coeff,"real");
    
    %solve the orthogonality equations to get the unknown coeff of ortho polys
    poly_coeff_sol=solve(ortho_eqn);
    
    if(N>2)  % when N=2 we have just one unkown and hence somehow the solve command returns a number instead of a structure (containing multiple numbers)
        poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';
    
    end
    %something to make output of solve command more usable (the output is a MATLAB structure by default, but we want to work with MATLAB sym array (symbolic array))
    %poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';
    
    %Substitute the values of unknown coeff in the ortho polys expressions
    polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);
    
    SP_polys_to_use=polys_to_use;
    
    %see the final ortho polys to be used 
    SP_polys_to_use';
    % Nc=N-1;
    thiag_polys_to_use=SP_polys_to_use(2:end);
    
    %see the final ortho polys to be used 
    % thiag_polys_to_use'
    rootP=vpasolve(thiag_polys_to_use(end));
    rootP=sort(abs(rootP(1:end/2)));
    rootP=[rootP' 1];
    
    %weights in the lagrange polynomial
    lgr_wts=sym("phi",[1 Nc+1]);
    Lgr_wts= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), lgr_wts, 'UniformOutput', false));
    
    for ii=1:length(rootP)
        Lgr_wts(ii)=1;
        for jj=1:length(rootP)
            if ~(ii==jj)
                Lgr_wts(ii)=Lgr_wts(ii)*(r^2-rootP(jj)^2)/(rootP(ii)^2-rootP(jj)^2);
            end
        end
    end

    quad_pts=rootP;
    quad_wts=sym("0");
    for ii=1:Nc+1
        %the 3 multiplier is taken into weights (from the integrand) as
        %done in thiagrajan et. al.
        quad_wts(ii)=int(Lgr_wts(ii)*3*r^2,r,0,1);
    end
end
% See the Model Equations

%use vpa (variable precision arithematic) to shorten and simplify the
%symbolic fractions often appearing in the galerkin Model
digits(dd)
Galerk_eqn_JP=vpa(Galerk_eqn_JP);
Galerk_eqn_SP=vpa(Galerk_eqn_SP);
Galerk_eqn_thiag=vpa(Galerk_eqn_thiag);
Galerk_eqn_thiag_GJ=vpa(Galerk_eqn_thiag_GJ);


%See the Galerkin Model
Galerk_eqn_JP;
Galerk_eqn_SP;
Galerk_eqn_thiag;
SM_2eqn;%%
dcdt=sym("dcdt",[1 Nc+1]);
Galerk_eqn_thiag_Explicit=Galerk_eqn_thiag;
for ii=1:Nc+1
    Galerk_eqn_thiag_Explicit=subs(Galerk_eqn_thiag_Explicit,diff(c(ii),t),dcdt(ii));
end
% digits(dd)
% Galerk_eqn_thiag_Explicit=vpasolve(Galerk_eqn_thiag_Explicit,dcdt);%% 
% Galerk_eqn_thiag_Explicit; Galerk_eqn_thiag_GJ; Test the Quadrature

%expression to be integrated
% expr_to_test=JP_polys_to_use(end)*SP_polys_to_use(end)*sin(r^2);

expr_to_test=expand(Lgr_wts(end)*diff(Lgr_wts(end),r));
% Test the GJ Quadrature

%see the points and weights being used in the quad
nn=10;
[pts,wts]=giveGaussJacobiQuad(nn,a,b);
int_quad_result=double(giveQuadResult(expr_to_test,nn,a,b));
int_exact_result=double(int(expr_to_test*r^2,r,0,1));
% Test the Thiag Quadrature

[pts_thiag,wts_thiag]=giveThiagQuad(nn);
int_quad_result=double(giveThiagQuadResult(expr_to_test,nn)/3);
int_exact_result=double(int(expr_to_test*r^2,r,0,1));

time_part_2=toc;
% Solve the Model Equations
% Solve the New Models
tic
%making the derived model that is in symbolic variables, to be usable by
%MATLABs ODE solvers
%M is the mass matrix, F is the RHS of the model, [M]*[d(alpha)dt]=[F]
[M_JP,F_JP]=massMatrixForm(Galerk_eqn_JP,alpha);
FF_JP=odeFunction(F_JP,alpha,delta(t));
M_JP=odeFunction(M_JP,alpha);

opt_JP=odeset("Mass",M_JP,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

[M_SP,F_SP]=massMatrixForm(Galerk_eqn_SP,alpha);
% F_SP
% M_SP
FF_SP=odeFunction(F_SP,alpha,delta(t));
M_SP=odeFunction(M_SP,alpha);

opt_SP=odeset("Mass",M_SP,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);


[M_thiag,F_thiag]=massMatrixForm(Galerk_eqn_thiag,c);
FF_thiag=odeFunction(F_thiag,c,delta(t));
M_thiag=odeFunction(M_thiag,c);

opt_thiag=odeset("Mass",M_thiag,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);


[M_thiag_GJ,F_thiag_GJ]=massMatrixForm(Galerk_eqn_thiag_GJ,c);
FF_thiag_GJ=odeFunction(F_thiag_GJ,c,delta(t));
M_thiag_GJ=odeFunction(M_thiag_GJ,c);

opt_thiag_GJ=odeset("Mass",M_thiag_GJ,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);


[M_SM,F_SM]=massMatrixForm(SM_2eqn,[c(1) c(2)]);
F_SM;
M_SM;
syms deltaDiff1(t) deltaDiff2(t) deltaDiff3(t)
F_SM=subs(F_SM,[diff(delta,t,3) diff(delta,t,2) diff(delta,t)],[deltaDiff3 deltaDiff2 deltaDiff1]);
% FF_SM=odeFunction(F_SM,c,delta(t));
FF_SM=odeFunction(F_SM,c,delta(t),deltaDiff1(t),deltaDiff2(t),deltaDiff3(t));
M_SM=odeFunction(M_SM,c);

opt_SM=odeset("Mass",M_SM,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);


% Calculate the initial values (at t=0) of alpha, such that conc. inside
% particle is uniform with value = per pStruct.C0 
for ii=1:N
    alpha0_JP(ii)=double(int(pStruct.C0*JP_polys_to_use(ii)*r^2,r,0,1)/int(JP_polys_to_use(ii)*JP_polys_to_use(ii)*r^2,r,0,1));
    alpha0_SP(ii)=double(int(pStruct.C0*SP_polys_to_use(ii)*r^2,r,0,1)/int(SP_polys_to_use(ii)*SP_polys_to_use(ii)*r^2,r,0,1));
    c0_thiag(ii)=pStruct.C0;
    c0_thiag_GJ(ii)=pStruct.C0;
end

c0_SM(1)=pStruct.C0;
c0_SM(2)=pStruct.C0;

%defining the dynamic flux 
flux = @(t) giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);

% [del,del1,del2,del3]=giveDeltafortau(5,pStruct.Deltavstauprofile,pStruct)
% [~,~,del2,~]=giveDeltafortau(5,pStruct.Deltavstauprofile,pStruct);


fluxVal = @(t) giveDeltaVal(t,pStruct);

fluxDiff1 = @(t) giveDeltaDiff1(t,pStruct);
fluxDiff2 = @(t) giveDeltaDiff2(t,pStruct);
fluxDiff3 = @(t) giveDeltaDiff3(t,pStruct);


function out = giveDeltaVal(t,pStruct)
    [out,~,~,~] = giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);
end
function out = giveDeltaDiff1(t,pStruct)
    [~,out,~,~] = giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);
end
function out = giveDeltaDiff2(t,pStruct)
    [~,~,out,~] = giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);
end
function out = giveDeltaDiff3(t,pStruct)
    [~,~,~,out] = giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);
end

time_part_3=toc;

[~,alpha_sol_JP]=ode15s(@(t,Alpha) FF_JP(t,Alpha,flux(t)),pStruct.taupoints,alpha0_JP,opt_JP);

SP_timer=tic;
[~,alpha_sol_SP]=ode15s(@(t,Alpha) FF_SP(t,Alpha,flux(t)),pStruct.taupoints,alpha0_SP,opt_SP);
SP_time=toc(SP_timer);

thiag_timer=tic;
[~,c_sol_thiag]=ode15s(@(t,CC) FF_thiag(t,CC,flux(t)),pStruct.taupoints,c0_thiag,opt_thiag);
thiag_time=toc(thiag_timer);

thiag_GJ_timer=tic;
[~,c_sol_thiag_GJ]=ode15s(@(t,CC) FF_thiag_GJ(t,CC,flux(t)),pStruct.taupoints,c0_thiag_GJ,opt_thiag_GJ);
thiag_GJ_time=toc(thiag_GJ_timer);

SM_timer=tic;
[~,c_sol_SM]=ode15s(@(t,CC) FF_SM(t,CC,fluxVal(t),fluxDiff1(t),fluxDiff2(t),fluxDiff3(t)),pStruct.taupoints,c0_SM,opt_SM);
SM_time=toc(SM_timer);

% Solve the PDEPE (benchmark) and Thiagrajan et.al. Galerkin Models
PDEPE_timer=tic;
[C_PDEPE, t_PDEPE] = givePDEPESol(rpoints,taupoints,pStruct);
PDEPE_time=toc(PDEPE_timer);

%% 
% clc

% figure;
label_fontsize = 8;
tick_fontsize = 12;
legend_fontsize = 12;
title_fontsize = 30; % New: font size for plot titles
fontname = 'Cambria Math';
interpreter = 'latex';
fig_width_cm = 16; % width in cm
fig_height_cm = 20; % height in cm (adjust as needed)
fig_width_in = fig_width_cm / 2.54; % convert to inches
fig_height_in = fig_height_cm / 2.54;

h_gap=0.13;
v_gap=0.1;
lf_mrg=0.13;
rt_mrg=0.02;
bt_mrg=0.07;
tp_mrg=0.05;

f1 = figure();
set(f1, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
set(f1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 fig_width_in fig_height_in]);
set(f1, 'PaperSize', [fig_width_in fig_height_in]);

set(gca, 'LooseInset', get(gca,'TightInset'))

Cs_sol_JP=zeros(1,length(taupoints))';
Cs_sol_SP=zeros(1,length(taupoints))';
SP_s = matlabFunction(subs(subs(Capprox_SP,r,1),alpha,aa));

for tt=1:length(taupoints)
    alpha_for_SP=num2cell(alpha_sol_SP(tt,:));
    Cs_sol_SP(tt)=SP_s(alpha_for_SP{:});
end

% plot(taupoints,Cs_sol_SP,"--m",'LineWidth', 2)
% set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
% xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
% ylabel('$C_s$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
% th = title('$\mathrm{Surface~Conc.~vs~Time}$', 'Interpreter', interpreter); % Set title handle
% set(th, 'FontSize', title_fontsize, 'FontName', fontname); % Explicitly set font size and font name

Cavg_sol_JP=zeros(1,length(taupoints));
Cavg_sol_SP=zeros(1,length(taupoints));
Cavg_thiag_new=zeros(1,length(taupoints));
Cavg_thiag_GJ=zeros(1,length(taupoints));

%create a matlab function out of the symbolic experssion (for faster computation)
JP_s = matlabFunction(subs(subs(Capprox_JP,r,1),alpha,aa));
SP_s = matlabFunction(subs(subs(Capprox_SP,r,1),alpha,aa));

Cavg_thiag_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));
Cavg_thiag_GJ_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));

%for thiag its just the last c, cs=c(end)
Cs_thiag_new=c_sol_thiag(:,end);
Cs_thiag_GJ=c_sol_thiag_GJ(:,end);

%for SM the variables themselves are just Cavg and Cs
Cavg_SM=c_sol_SM(:,1);
Cs_SM=c_sol_SM(:,2);

for tt=1:length(taupoints)
    %surface conc.
    alpha_for_JP=num2cell(alpha_sol_JP(tt,:));
    Cs_sol_JP(tt)=JP_s(alpha_for_JP{:});
    alpha_for_SP=num2cell(alpha_sol_SP(tt,:));
    Cs_sol_SP(tt)=SP_s(alpha_for_SP{:});

    % average conc. equal to alpha_1
    Cavg_sol_JP(tt)=alpha_sol_JP(tt,1);
    Cavg_sol_SP(tt)=alpha_sol_SP(tt,1);

    c_for_thiag=num2cell(c_sol_thiag(tt,:));
    Cavg_thiag_new(tt)=Cavg_thiag_fn(c_for_thiag{:});

    c_for_thiag_GJ=num2cell(c_sol_thiag_GJ(tt,:));
    Cavg_thiag_GJ(tt)=Cavg_thiag_GJ_fn(c_for_thiag_GJ{:});
end
% Calculate the average conc. for PDEPE solution

Cavg_PDEPE=zeros(1,length(taupoints));
for ii=1:length(taupoints)
    Cr2=C_PDEPE(ii,:).*rpoints.^2;
    Cavg_PDEPE(ii)=(3)*trapz(rpoints,Cr2);
end

CC=linspace(0,1,100);
f_at_C=zeros(1,length(CC));
for ii=1:length(CC)
    f_at_C(ii)=giveFforC(CC(ii),pStruct.FvsCprofile,pStruct);
end
%  Plot the difference between Models

model_names = {...
    ['$\mathbf{PDEPE}-$', num2str(length(pStruct.rpoints)), '$\;\mathrm{nodes}$'], ...
    ['$\mathbf{$', num2str(Nc+1), '$}\;\mathrm{Eqn}-\mathbf{ThiagOriginal}-$', num2str(Nc), '$\;\mathrm{QdPts}$'], ...
    ['$\mathbf{$', num2str(Nc+1), '$}\;\mathrm{Eqn}-\mathbf{ThiagModified}-$', num2str(N_quad_thiag_GJ), '$\;\mathrm{QdPts}$'], ...
    ['$\mathbf{$', num2str(N), '$}\;\mathrm{Eqn}-\mathbf{New}-$', num2str(N_quad_SP), '$\;\mathrm{QdPts}$'] ...
    };

Cs_All=[C_PDEPE(:,end), Cs_thiag_new, Cs_thiag_GJ, Cs_sol_SP];
Cavg_All = [Cavg_PDEPE; Cavg_thiag_new; Cavg_thiag_GJ; Cavg_sol_SP];

Cavg_Analytic=zeros(1,length(taupoints));
% Cavg_Analytic(1)=pStruct.C0;
% 
flux_at_tau = zeros(1, length(taupoints));
for ii = 1:length(taupoints)
    flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
end
% 
% 
% for ii=2:length(taupoints)
%     Cavg_Analytic(ii) = pStruct.C0-(3)*trapz(taupoints(1:ii),flux_at_tau(1:ii));
% end

taupoints_for_Cavg_calc=linspace(taupoints(1),taupoints(end),20*length(taupoints));
Cavg_Analytic_dummy=zeros(1,length(taupoints_for_Cavg_calc));
Cavg_Analytic_dummy(1)=pStruct.C0;



% for ii=2:length(taupoints)
%     Cavg_Analytic(ii) = pStruct.C0-(3)*trapz(taupoints(1:ii),flux_at_tau(1:ii));
% end

flux_at_tau_for_Cavg_calc = zeros(1, length(taupoints_for_Cavg_calc));
for ii = 1:length(taupoints_for_Cavg_calc)
    flux_at_tau_for_Cavg_calc(ii) = giveDeltafortau(taupoints_for_Cavg_calc(ii), pStruct.Deltavstauprofile, pStruct);
end

for ii=2:length(taupoints_for_Cavg_calc)
    Cavg_Analytic_dummy(ii) = pStruct.C0-(3)*trapz(taupoints_for_Cavg_calc(1:ii),flux_at_tau_for_Cavg_calc(1:ii));
end

Cavg_Analytic=interp1(taupoints_for_Cavg_calc,Cavg_Analytic_dummy,taupoints);


clr_used=["blue", "black", "green", "magenta"];

[~, cols]=size(Cs_All);
num_models=cols;

taus_idx = [ceil((1/5) * length(taupoints)), 2*377, ...
            2*544, ceil((4/5) * length(taupoints))];

profile_plt_loc = [15 16 21 22];
lnWdt=1.5;

num_pairs = nchoosek(num_models, 2);
pair_idx = 1;
iplot = 7;

model_names_old=["\bfPDEPE\rm-\bf"+length(pStruct.rpoints)+"\rmnodes",...
    "\bf"+(Nc+1)+"\rmEqn"+"-\bfThiagOriginal\rm-\bf"+Nc+"\rmQdPts",...
    "\bf"+(Nc+1)+"\rmEqn"+"-\bfThiagModified\rm-\bf"+N_quad_thiag_GJ+"\rmQdPts",...
    "\bf"+N+"\rmEqn"+"-\bfNewModal\rm-\bf"+N_quad_SP+"\rmQdPts"];

Cs_err_with_time=zeros(6,length(taupoints));
Cavg_err_with_time=zeros(6,length(taupoints));

for kk=[9,10]
    for ii = 1:num_models
        for jj = ii+1:num_models
 
            Cs_err_with_time(pair_idx,:)=Cs_All(:,ii)-Cs_All(:,jj);
            % Cavg_err_with_time(pair_idx,:)=Cavg_All(ii,:)-Cavg_All(jj,:);
            % Cavg_err_with_time(pair_idx,:)=Cavg_Analytic-Cavg_All(jj,:);

    
            pair_idx=pair_idx+1;
        end
    end
    
    for ii = 1:num_models
            Cavg_err_with_time(ii,:)=Cavg_Analytic-Cavg_All(ii,:);
    end

    subtightplot(5, 2, kk, [v_gap h_gap], [bt_mrg tp_mrg], [lf_mrg rt_mrg]);
    if(kk==9)
        plot(taupoints,Cs_err_with_time(1,:),'-k', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cs_err_with_time(2,:),'-g', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cs_err_with_time(3,:),'--m', 'LineWidth', 1 ,'Marker', 'none');
        xlabel('$\tau$', 'FontSize', label_fontsize-4, 'FontName', fontname, 'Interpreter', interpreter);
        ylabel('${C_{\mathrm{s,PDEPE}}}-{C_{\mathrm{s,Model}}}$', 'FontSize', label_fontsize-4, 'FontName', fontname, ...
            'Interpreter', interpreter);
        th = title("Error in Surface Conc."); % Set title handle
        set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
        set(gca, 'FontName', fontname, 'FontSize', tick_fontsize, 'TickLabelInterpreter', interpreter);
    else
        % plot(taupoints,Cavg_err_with_time(1,:),'-b', 'LineWidth', 1 ,'Marker', 'none');
        % hold on;
        % plot(taupoints,Cavg_err_with_time(2,:),'-k', 'LineWidth', 1 ,'Marker', 'none');
        % hold on;
        % plot(taupoints,Cavg_err_with_time(3,:),'-g', 'LineWidth', 1 ,'Marker', 'none');
        % hold on;
        % plot(taupoints,Cavg_err_with_time(4,:),'--m', 'LineWidth', 1 ,'Marker', 'none');
        % hold on;

        semilogy(taupoints,abs(Cavg_err_with_time(1,:))+1e-10,'-b', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        semilogy(taupoints,abs(Cavg_err_with_time(2,:))+1e-10,'-k', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        semilogy(taupoints,abs(Cavg_err_with_time(3,:))+1e-10,'-g', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        semilogy(taupoints,abs(Cavg_err_with_time(4,:))+1e-10,'--m', 'LineWidth', 1 ,'Marker', 'none');
        hold on;


        xlabel('$\tau$', 'FontSize', label_fontsize-4, 'FontName', fontname, 'Interpreter', interpreter);
        ylabel('$\left|{C_{\mathrm{avg,Anly}}}-{C_{\mathrm{avg,Model}}}\right|+10^{-10}$', 'FontSize', label_fontsize-4, 'FontName', fontname, ...
            'Interpreter', interpreter);
        th = title("Error in Average Conc."); % Set title handle
        set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
        set(gca, 'FontName', fontname, 'FontSize', 10, 'TickLabelInterpreter', interpreter);
        % ax=gca;
        yticks([1e-12 1e-10 1e-8 1e-6 1e-4])
        yticklabels({"$10^{-12}$","$10^{-10}$","$10^{-8}$","$10^{-6}$"})
        ax = gca;
        ax.YAxis.Exponent = -7;

        % text(ax.Position(1)-0.05, ax.Position(2)+ax.Position(4), '×10^{-7}', ...
        % 'FontSize', 30, 'Interpreter', 'tex', 'HorizontalAlignment', 'left');

        % ax=gca;
        % ax.YAxis.ExponentLocation = 'left';

    end
end

subtightplot(5,2,[3 4 5 6],[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])

yyaxis right

% Cavg_Analytic=zeros(1,length(taupoints));
% Cavg_Analytic(1)=pStruct.C0;
% for ii=2:length(taupoints)
%     Cavg_Analytic(ii) = pStruct.C0-(3)*trapz(taupoints(1:ii),flux_at_tau(1:ii));
% end
% figure
% plot(taupoints,Cavg_Analytic)


plot(taupoints, flux_at_tau, ':k', 'LineWidth', 1.1, 'Color', [0 0 0 0.8])
yticks([])
ax = gca;
ax.YColor = 'k';

yyaxis left
plot(taupoints,C_PDEPE(:,end),"b",'LineWidth', lnWdt)
hold on
plot(taupoints,Cs_thiag_new,"-k",'LineWidth', lnWdt-0.5)
hold on
plot(taupoints,Cs_thiag_GJ,"-g",'LineWidth', lnWdt-0.5)
hold on
plot(taupoints,Cs_sol_SP,"--m",'LineWidth', lnWdt-0.5)

xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$C_s$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Surface~Conc.~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
legend_handle = legend(model_names_old, 'Location', 'northwest', 'FontSize', legend_fontsize, ...
    'FontName', fontname);
ax = gca;
ax.YColor = 'k';
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
ylim([-0.1,1])

subtightplot(5,2,8,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
model_runtimes = [PDEPE_time, thiag_time, thiag_GJ_time, SP_time];
bar_handle = bar(model_runtimes, 0.5,'FaceColor', 'flat');
set(gca, 'FontSize', tick_fontsize, 'FontName', fontname);

bar_handle.CData = [
    0, 0, 1;       % blue
    0, 0, 0;       % black
    0, 1, 0;       % green
    1, 0, 1;       % magenta
];
set(gca,'YScale', 'log');
ylabel('$\mathrm{Runtime~[s]}$', 'FontSize', label_fontsize-6, 'FontName', fontname, 'Interpreter', interpreter);
th = title('$\mathrm{Model~Runtimes}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
ylim([1e-3 1e3]);
for ii = 1:length(model_runtimes)
    text(ii, model_runtimes(ii)*1.1, sprintf('$%.0f\\mathrm{~ms}$', model_runtimes(ii) * 1000), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, ...
        'FontName', fontname, 'Interpreter', interpreter);
end
set(gca,'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'FontName', fontname, 'FontSize', tick_fontsize);
xticklabels([])

yticks([1e-2 1 1e2])
yticklabels({'$10^{-2}$', '$1$', '$10^{2}$'})
set(gca, 'TickLabelInterpreter', interpreter, 'FontSize', tick_fontsize, 'FontName', fontname);


subtightplot(5,2,7,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
plot(taupoints,Cavg_Analytic,"--r")
hold on
plot(taupoints,Cavg_PDEPE,"b")
hold on
plot(taupoints,Cavg_thiag_new,"--k")
hold on
plot(taupoints,Cavg_thiag_GJ,"-.g")
hold on
plot(taupoints,Cavg_sol_SP,"--m")

xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$C_{\mathrm{avg}}$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Average~Conc.~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)

subtightplot(5,2,2,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
% flux_at_tau = zeros(1, length(taupoints));
% for ii = 1:length(taupoints)
%     flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
% end
plot(taupoints,flux_at_tau)
xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$\delta$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Flux~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
ylim([-0.04 0.08])
yticks([-0.04 0 0.04])
yticklabels({'-0.04', '0', '0.04'})

subtightplot(5,2,1,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
semilogy(CC,f_at_C)
xlabel('$C$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$f$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Diffusivity~vs~Conc.}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
yticks([1e-3 1e-2 1e-1])
yticklabels({'$10^{-3}$', '$10^{-2}$', '$10^{-1}$'})
set(gca, 'TickLabelInterpreter', interpreter, 'FontSize', tick_fontsize, 'FontName', fontname);

% flux_at_tau=zeros(1,length(taupoints));
% for ii=1:length(taupoints)
%     flux_at_tau(ii)=giveDeltafortau(taupoints(ii),pStruct.Deltavstauprofile,pStruct);
% end
%% 


label_fontsize = 8;
tick_fontsize = 12;
legend_fontsize = 12;
title_fontsize = 30; % New: font size for plot titles
fontname = 'Cambria Math';
interpreter = 'latex';
fig_width_cm = 16; % width in cm
fig_height_cm = 20; % height in cm (adjust as needed)
fig_width_in = fig_width_cm / 2.54; % convert to inches
fig_height_in = fig_height_cm / 2.54;

h_gap=0.13;
v_gap=0.1;
lf_mrg=0.13;
rt_mrg=0.02;
bt_mrg=0.07;
tp_mrg=0.05;

f2=figure();

set(f2, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
set(f2, 'PaperUnits', 'inches', 'PaperPosition', [0 0 fig_width_in fig_height_in]);
set(f2, 'PaperSize', [fig_width_in fig_height_in]);

set(gca, 'LooseInset', get(gca,'TightInset'))

taus_idx = [ceil((1/5) * length(taupoints)), 2*377, ...
            2*544, ceil((4/5) * length(taupoints))];

subtightplot(3,2,[1 2],[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
plot(taupoints,flux_at_tau)
xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$\delta$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Flux~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
ylim([-0.04 0.08])

for ii=1:length(taus_idx)
    hold on
    xline(taupoints(taus_idx(ii)),"--",LineWidth=0.5)
end


% Plot the inside Conc.Profile in the specified subplot locations
profile_plt_loc = [3 4 5 6];
lnWdt=1.5;
for idx = 1:4
    subtightplot(3, 2, profile_plt_loc(idx),[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
    plot(rpoints, C_PDEPE(taus_idx(idx), :), "b",'LineWidth', lnWdt);
    hold on;
    C_thiag_new_vs_r = double(subs(subs(Capprox_thiag, c, c_sol_thiag(taus_idx(idx), :)), r, rpoints));
    plot(rpoints, C_thiag_new_vs_r, "-k",'LineWidth', lnWdt-0.5);
    hold on;
    C_thiag_GJ_vs_r = double(subs(subs(Capprox_thiag, c, c_sol_thiag_GJ(taus_idx(idx), :)), r, rpoints));
    plot(rpoints, C_thiag_GJ_vs_r, "-g",'LineWidth', lnWdt-0.5);
    hold on;
    C_SP_vs_r = double(subs(subs(Capprox_SP, alpha, alpha_sol_SP(taus_idx(idx), :)), r, rpoints));
    plot(rpoints, C_SP_vs_r, "--m",'LineWidth', lnWdt-0.5);
    ylim([0 1.5]);
    xlabel("\it{r}", 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
    ylabel("\it{C}", 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
    % title("\it{C} vs \it{r} ($\tau= $" + taupoints(taus_idx(idx))+")",'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
    title("Internal Conc. ($\tau =$ " + sprintf("%.1f",taupoints(taus_idx(idx)))+")",'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);

    set(gca, 'TickLabelInterpreter', interpreter, 'FontSize', tick_fontsize, 'FontName', fontname);
    
    if(idx==4)
        legend_handle = legend(model_names_old, 'Location', 'northwest', 'FontSize', legend_fontsize-3, ...
    'FontName', fontname);
    end
end