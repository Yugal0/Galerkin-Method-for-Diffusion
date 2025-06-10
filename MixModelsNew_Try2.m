%% Creating Variables
clear all
clc
dd=8
digits(dd)

N=6;
N_quad=14;

function [Capprox_thiag,Lgr_wts]=giveCapproxThiag(N)
    syms r t
    assume(r,"real")
    assume(t,"real")

    p=sym("P",[1 N]);
    P= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p, 'UniformOutput', false));

    aa=sym("alpha",[1 N]);
    alpha= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa, 'UniformOutput', false));

    Capprox = sym("0");  
    for ii=1:N    
        Capprox=Capprox+P(ii)*alpha(ii);
    end
    
    polys_to_use=sym("0");

    poly_coeff=sym("nu",[N N]);

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
    polys_to_use';

    kk=1;
    poly_coeff_in_use=sym("0");
    for ii=2:N
    for jj=1:ii-1
        poly_coeff_in_use(kk)=poly_coeff(ii,jj);
        kk=kk+1;
    end
    end

    ortho_eqn=sym("0");
    kk=1;
    for ii=1:N
    for jj=1:N
        if(ii~=jj && ii<jj)           
            ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
            kk=kk+1;
        end
    end
    end

    assume(poly_coeff,"real");

    poly_coeff_sol=solve(ortho_eqn);

    if(N>2)  
        poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';
    end

    polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);

    SP_polys_to_use=polys_to_use;

    SP_polys_to_use';

    Capprox_SP=subs(Capprox,P,SP_polys_to_use);

    Coeff_SP_poly = coeffs(Capprox_SP,r);
    Coeff_SP_poly = subs(Coeff_SP_poly,alpha,aa);
    Nc=N-1;
    thiag_polys_to_use=SP_polys_to_use(2:end);
    
    thiag_polys_to_use';
    rootP=vpasolve(thiag_polys_to_use(end));
    rootP=sort(abs(rootP(1:end/2)));
    rootP=[rootP' 1];
    
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
    cc = sym("c",[1 Nc+1]);
    c = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), cc, 'UniformOutput', false));
    
    Capprox_thiag=subs(Capprox,[P alpha],[Lgr_wts c]);
end

function [Capprox_SP,SP_polys_to_use]=giveCapproxNew(N)
    syms r t
    assume(r,"real")
    assume(t,"real")

    p=sym("P",[1 N]);
    P= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p, 'UniformOutput', false));

    aa=sym("alpha",[1 N]);
    alpha= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa, 'UniformOutput', false));

    Capprox = sym("0");  
    for ii=1:N    
    Capprox=Capprox+P(ii)*alpha(ii);
    end

    polys_to_use=sym("0");

    poly_coeff=sym("nu",[N N]);

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
    polys_to_use';

    kk=1;
    poly_coeff_in_use=sym("0");
    for ii=2:N
    for jj=1:ii-1
        poly_coeff_in_use(kk)=poly_coeff(ii,jj);
        kk=kk+1;
    end
    end

    ortho_eqn=sym("0");
    kk=1;
    for ii=1:N
    for jj=1:N
        if(ii~=jj && ii<jj)           
            ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
            kk=kk+1;
        end
    end
    end

    assume(poly_coeff,"real");

    poly_coeff_sol=solve(ortho_eqn);

    if(N>2)  
    poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';

    end

    polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);

    SP_polys_to_use=polys_to_use;

    SP_polys_to_use';

    Capprox_SP=subs(Capprox,P,SP_polys_to_use);

    Coeff_SP_poly = coeffs(Capprox_SP,r);
    Coeff_SP_poly = subs(Coeff_SP_poly,alpha,aa);
end

function [Galerk_eqn_thiag]=giveGalerkinThiagOrg(N,pStruct)
    syms r t delta(t)
    assume(r,"real")
    assume(t,"real")
    [Capprox_thiag, Lgr_wts]=giveCapproxThiag(N);
    W_thiag=Lgr_wts;
    N_quad_thiag=N-1;
    ff = @(r,Capprox) firstOutputOnly(r,Capprox,pStruct);
    for ii=1:N
        Galerk_eqn_thiag(ii)=int(W_thiag(ii)*3*r^2*diff(Capprox_thiag,t),r,0,1)==subs(W_thiag(ii)*3*r^2*(-delta),r,1)-giveThiagQuadResult(diff(W_thiag(ii),r)*ff(r,Capprox_thiag)*diff(Capprox_thiag,r),N_quad_thiag);
    end
end

function [Galerk_eqn_thiag_GJ]=giveGalerkinThiagModify(N,N_quad_thiag_GJ,pStruct)
    syms r t delta(t)
    assume(r,"real")
    assume(t,"real")
    [Capprox_thiag, Lgr_wts]=giveCapproxThiag(N);
    W_thiag=Lgr_wts;
    N_quad_thiag=N-1;
    ff = @(r,Capprox) firstOutputOnly(r,Capprox,pStruct);
    a=0;
    b=2;
    for ii=1:N
        Galerk_eqn_thiag_GJ(ii)=int(W_thiag(ii)*3*r^2*diff(Capprox_thiag,t),r,0,1)==subs(W_thiag(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_thiag(ii),r)*3*ff(r,Capprox_thiag)*diff(Capprox_thiag,r),N_quad_thiag_GJ,a,b);
    end
end

function [Galerk_eqn_SP]=giveGalerkinNew(N,N_quad_SP,pStruct)
    syms r t delta(t)
    assume(r,"real")
    assume(t,"real")
    [Capprox_SP, SP_polys_to_use]=giveCapproxNew(N);
    W_SP=SP_polys_to_use;
    % N_quad_thiag=N-1;   
    ff = @(r,Capprox) firstOutputOnly(r,Capprox,pStruct);
    a=0;
    b=2;   
    for ii=1:N
        Galerk_eqn_SP(ii)=int(W_SP(ii)*3*r^2*diff(Capprox_SP,t),r,0,1)==subs(W_SP(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_SP(ii),r)*3*ff(r,Capprox_SP)*diff(Capprox_SP,r),N_quad_SP,a,b);
    end
end

pStruct.FvsCprofile=char("nmc811AMIDRFitted5");
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


N_quad_thiag=14;
Galerk_eqn_thiag=giveGalerkinThiagOrg(N,pStruct);
Galerk_eqn_thiag_GJ=giveGalerkinThiagModify(N,N_quad_thiag,pStruct);
Galerk_eqn_SP=giveGalerkinNew(N,N_quad,pStruct);

N_low=3;
N_high=6;
N_quad_low=8;
N_quad_high=14;
Galerk_eqn_SP_low=giveGalerkinNew(N_low,N_quad_low,pStruct);
Galerk_eqn_SP_high=giveGalerkinNew(N_high,N_quad_high,pStruct);

% digits(dd)
Galerk_eqn_SP=vpa(Galerk_eqn_SP);
Galerk_eqn_thiag=vpa(Galerk_eqn_thiag);
Galerk_eqn_thiag_GJ=vpa(Galerk_eqn_thiag_GJ);
Galerk_eqn_SP_low=vpa(Galerk_eqn_SP_low);
Galerk_eqn_SP_high=vpa(Galerk_eqn_SP_high);

syms r t delta(t)
assume(r,"real")
assume(t,"real")
p = sym("P",[1 N]);
P = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p, 'UniformOutput', false));
aa = sym("alpha",[1 N]);
alpha = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa, 'UniformOutput', false));
cc=sym("c",[1 N]);
c= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), cc, 'UniformOutput', false));

[M_thiag,F_thiag]=massMatrixForm(Galerk_eqn_thiag,c);
FF_thiag=odeFunction(F_thiag,c,delta(t));
M_thiag=odeFunction(M_thiag,c);
opt_thiag=odeset("Mass",M_thiag,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

[M_thiag_GJ,F_thiag_GJ]=massMatrixForm(Galerk_eqn_thiag_GJ,c);
FF_thiag_GJ=odeFunction(F_thiag_GJ,c,delta(t));
M_thiag_GJ=odeFunction(M_thiag_GJ,c);
opt_thiag_GJ=odeset("Mass",M_thiag_GJ,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

[M_SP,F_SP]=massMatrixForm(Galerk_eqn_SP,alpha);
FF_SP=odeFunction(F_SP,alpha,delta(t));
MM_SP=odeFunction(M_SP,alpha);
opt_SP=odeset("Mass",MM_SP,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

p_low = sym("P",[1 N_low]);
P_low = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p_low, 'UniformOutput', false));
aa_low = sym("alpha",[1 N_low]);
alpha_low = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa_low, 'UniformOutput', false));

[M_SP_low,F_SP_low]=massMatrixForm(Galerk_eqn_SP_low,alpha_low);
FF_SP_low=odeFunction(F_SP_low,alpha_low,delta(t));
MM_SP_low=odeFunction(M_SP_low,alpha_low);
opt_SP_low=odeset("Mass",MM_SP_low,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

p_high = sym("P",[1 N_high]);
P_high = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(r)'), r), p_high, 'UniformOutput', false));
aa_high = sym("alpha",[1 N_high]);
alpha_high = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa_high, 'UniformOutput', false));

[M_SP_high,F_SP_high]=massMatrixForm(Galerk_eqn_SP_high,alpha_high);
FF_SP_high=odeFunction(F_SP_high,alpha_high,delta(t));
MM_SP_high=odeFunction(M_SP_high,alpha_high);
opt_SP_high=odeset("Mass",MM_SP_high,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

%% Solving the models
%set inital values (acc projection method)
[~, SP_polys_to_use]=giveCapproxNew(N);
[~, SP_polys_to_use_low]=giveCapproxNew(N_low);
[Capprox_SP_high, SP_polys_to_use_high]=giveCapproxNew(N_high);
[Capprox_SP_low, SP_polys_to_use_low]=giveCapproxNew(N_low);

for ii=1:N
    alpha0_SP(ii)=double(int(pStruct.C0*SP_polys_to_use(ii)*r^2,r,0,1)/int(SP_polys_to_use(ii)*SP_polys_to_use(ii)*r^2,r,0,1));
    c0_thiag(ii)=pStruct.C0;
    c0_thiag_GJ(ii)=pStruct.C0;
end

for ii=1:N_low
    alpha0_SP_low(ii)=double(int(pStruct.C0*SP_polys_to_use_low(ii)*r^2,r,0,1)/int(SP_polys_to_use_low(ii)*SP_polys_to_use_low(ii)*r^2,r,0,1));
end
for ii=1:N_high
    alpha0_SP_high(ii)=double(int(pStruct.C0*SP_polys_to_use_high(ii)*r^2,r,0,1)/int(SP_polys_to_use_high(ii)*SP_polys_to_use_high(ii)*r^2,r,0,1));
end

%defining the dynamic flux 
flux = @(t) giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);

fluxVal = @(t) giveDeltaVal(t,pStruct);

fluxDiff1 = @(t) giveDeltaDiff1(t,pStruct);
fluxDiff2 = @(t) giveDeltaDiff2(t,pStruct);
fluxDiff3 = @(t) giveDeltaDiff3(t,pStruct);

thiag_timer=tic;
[~,c_sol_thiag]=ode15s(@(t,CC) FF_thiag(t,CC,flux(t)),pStruct.taupoints,c0_thiag,opt_thiag);
thiag_time=toc(thiag_timer);

thiag_GJ_timer=tic;
[~,c_sol_thiag_GJ]=ode15s(@(t,CC) FF_thiag_GJ(t,CC,flux(t)),pStruct.taupoints,c0_thiag_GJ,opt_thiag_GJ);
thiag_GJ_time=toc(thiag_GJ_timer);

SP_timer=tic;
[~,alpha_sol_SP]=ode15s(@(t,Alpha) FF_SP(t,Alpha,flux(t)),pStruct.taupoints,alpha0_SP,opt_SP);
SP_time=toc(SP_timer);

PDEPE_timer=tic;
[C_PDEPE, t_PDEPE] = givePDEPESol(rpoints,taupoints,pStruct);
% [C_PDEPE, t_PDEPE] = [zeros(length(taupoints),length(taupoints)) taupoints]
% C_PDEPE=zeros(length(taupoints),length(rpoints));
% t_PDEPE=taupoints;
PDEPE_time=toc(PDEPE_timer);

%%
clc

syms t delta(t)
aa_high = sym("alpha",[1 length(F_SP_high)]);
alpha_high = cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa_high, 'UniformOutput', false));

aa_low= sym("alpha",[1 length(F_SP_low)]);
alpha_low= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), aa_low, 'UniformOutput', false));

F_SP_low_new=sym(zeros(1,length(F_SP_high)));
F_SP_low_new(1,1:length(F_SP_low))=F_SP_low;
% F_SP_low_new(1,length(F_SP_low)+1:end)=alpha_high(length(F_SP_low)+1:end);

F_SP_low_new'

FF_SP_low=odeFunction(F_SP_low_new,alpha_high,delta(t));

FF_SP_low_2=odeFunction(F_SP_low,alpha_low,delta(t));

FF_SP_high=odeFunction(F_SP_high,alpha_high,delta(t));

M_SP_low_new = eye(length(M_SP_high));
M_SP_low_new(1:length(M_SP_low),1:length(M_SP_low))=M_SP_low;

M_SP_low_new

MM_SP_low=odeFunction(M_SP_low_new,alpha_low);

MM_SP_low_2=odeFunction(M_SP_low,alpha_low);

MM_SP_high=odeFunction(M_SP_high,alpha_high);

opt_SP_mix=odeset("Mass",@(t,alpha) conditionalMassMatrix(t,alpha,flux(t),MM_SP_low,MM_SP_high,N_low),"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);

Mix_timer=tic;
[~,alpha_sol_SP_mix]=ode15s(@(t,Alpha) FF_SP_mix(t,Alpha,flux(t),FF_SP_low,FF_SP_high,N_low),pStruct.taupoints,alpha0_SP_high,opt_SP_mix);
Mix_time=toc(Mix_timer);
%% 
ff_new = @(Capprox) firstOutputOnly(r,Capprox,pStruct);
[Capprox_low, poly_low]=giveCapproxNew(N_low);
[Capprox_high, poly_high]=giveCapproxNew(N_high);

ff_low=matlabFunction(ff_new(subs(Capprox_low,[alpha_low r],[aa_low 1])));
ff_high=matlabFunction(ff_new(subs(Capprox_high,[alpha_high r],[aa_high 1])));

Csurf_low_fun=matlabFunction(subs(Capprox_low,[alpha_low r],[aa_low 1]))
Csurf_high_fun=matlabFunction(subs(Capprox_high,[alpha_high r],[aa_high 1]))

poly_surf_low=subs(poly_low,[alpha_low r],[aa_low 1])
poly_surf_high=subs(poly_high,[alpha_low r],[aa_low 1])


%% 

clc
% 
opts_high=odeset("Events",@(t,alpha) event_high(t, alpha,flux(t),ff_high),"Mass",MM_SP_high(t,alpha),"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);
opts_low=odeset("Events",@(t,alpha) event_low(t, alpha,flux(t),ff_low),"Mass",MM_SP_low_2(t,alpha),"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);
opts_to_use=opts_low;

alpha0_to_use=alpha0_SP_low;

% alpha0_to_use=alpha0_SP_low;
% alpha0_to_use=[1 1 1 1 1 1];

% tspan=pStruct.taupoints;

% RHS_to_use=FF_SP_high;

RHS_to_use=FF_SP_low_2;
% RHS_to_use

tspan_full=pStruct.taupoints;
tf=tspan_full(end);
tspan=tspan_full;
te=0;

alpha_sol_total=[alpha0_to_use zeros(1,N_high-N_low)];
trigger_count=0;
tspan_run=[0];
te_list=[];
Mix_New_timer=tic;
while te<tf
    [t_run,alpha_sol_present,te,alpha_e,ie]=ode15s(@(t,Alpha) RHS_to_use(t,Alpha,flux(t)),tspan,alpha0_to_use,opts_to_use);
    trigger_count=trigger_count+ie;
    
    te_list=[te_list; te];
    tspan_run=[tspan_run; t_run(2:end)];

    if(isempty(te))
        te=tspan_full(end);
        % break;
    end
    
    % [~,idx_of_te]=min(abs(tspan_full-te));

    smaller_indices = find(tspan_full <= te);
    if ~isempty(smaller_indices)
        idx_of_te = smaller_indices(end); 
    else
        idx_of_te = NaN;  
    end

    if(length(alpha_sol_present(end,:))<N_high)
        dummy_loc=zeros(length(alpha_sol_present),N_high);
        dummy_loc(:,1:N_low)=alpha_sol_present;
        alpha_sol_present=dummy_loc;
    end
    alpha_sol_present
    size(alpha_sol_present)
    alpha_sol_total=[alpha_sol_total;alpha_sol_present(2:end,:)];
    
    if(te==tspan_full(end))
        break;
    end

    if(~mod(trigger_count,2))
        display("swicthed to N="+N_low+" at t="+te);
        RHS_to_use=FF_SP_low_2;
        opts_to_use=opts_low;
        alpha0_to_use=alpha_sol_present(end,1:N_low);
        alpha0_last=[zeros(1,N_low-1), alpha_sol_present(end,N_low:end)];
        alpha0_last=num2cell(alpha0_last);
        alpha0_to_use(end)=Csurf_high_fun(alpha0_last{:})/vpa(poly_surf_low(N_low));
    else
        display("swicthed to N="+N_high+" at t="+te);
        RHS_to_use=FF_SP_high;
        opts_to_use=opts_high;
        alpha0_to_use=alpha_sol_present(end,:);
        % alpha0_to_use
    end

    tspan=linspace(te,tspan_full(end),length(tspan_full)-idx_of_te);
    % tspan=tspan_full(idx_of_te:end);
end

Mix_New_time=toc(Mix_New_timer);

function [value, isterminal, direction] = event_high(t,alpha,flux,ff_high)
    % trigger_val = t-2;
    alpha=num2cell(alpha);
    trigger_val = (abs(-flux/ff_high(alpha{:}))-1);

    value = trigger_val;  
    isterminal = 1;                
    direction = -1;                
end

function [value, isterminal, direction] = event_low(t,alpha,flux,ff_low)
    % trigger_val = t-2;
    alpha=num2cell(alpha);
    trigger_val = (abs(-flux/ff_low(alpha{1:length(alpha)}))-1);
    value = trigger_val;  
    isterminal = 1;                
    direction = 1;                
end


%%%
% ode_fun_to_use=FF_SP_high;
% opt_SP_mix=odeset("Mass",MM_SP_high(t,alpha),"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);
% 
% while t0<tf 
%     tspan=[t0 tf];
%     options=odeset("Events",@my_events);
% 
%     [t,alpha_sol_mix,te,alpha_e,ie]=ode15s(@(t,Alpha) ode_fun_to_use(flux),pStruct.taupoints,alpha0_SP_high,opt_SP_mix);
% 
%     t_total = [t_total; t];
%     alpha_sol_total = [alpha_sol_total; alpha_sol_mix];
% 
%     if(ie==1)
%         alpha_0=[1, zeros(1,N_low)];
%         ode_fun_to_use=FF_SP_low_2;
%         mass_fun_to_use=;
%     else
%         ode_fun_to_use=FF_SP_high;
%         mass_fun_to_use=;
%     end
% 
%     t0=te+1e-6;
% end
% 
% function f=FF_Mix_wrapper(t,Alpha,flux)
%     if(trigger is there) 
%         f=@FF_SP_low;
%     else
%         f=@FF_SP_high;
%     end
% end

% function [value, isterminal, direction] = my_events(t, c)
%     trigger_val = t-2;
%     value = trigger_val;  % Two events, same condition
%     isterminal = 1;                 % Stop on either
%     direction = 1;                 % Event 1: upward, Event 2: downward
% end


function FF = FF_SP_mix(tau,alpha,flux,FF_SP_low,FF_SP_high,N_low)    
    % if(abs(alpha(3)/alpha(2))>0.1 && alpha(3)>0 && alpha(2)>0)
    % if(abs(flux)>0.03)
    if(tau>2)
        FF=FF_SP_high(tau,alpha,flux);
    else
        % alpha(N_low+1:end)=0;
        FF=FF_SP_low(tau,alpha,flux);
    end
end

function MM = conditionalMassMatrix(tau,alpha,flux,MM_SP_low,MM_SP_high,N_low)
    % if(abs(alpha(3)/alpha(2))>0.1 && alpha(3)>0 && alpha(2)>0)
    % if(abs(flux)>0.01)
    if(tau>2)
        MM=MM_SP_high(tau,alpha);
    else
        % assignin("base","alpha_vall",alpha)
        % alpha(N_low+1:end)=0;
        MM=MM_SP_low(tau,alpha);
    end
end

%% 
function out1 = firstOutputOnly(r,Capprox,pStruct)
    [out1, ~] = giveFforC(Capprox, pStruct.FvsCprofile, pStruct);
end

function result = giveQuadResult(integrand,n,a,b)
    syms r
    [quad_pts,quad_wts]=giveGaussJacobiQuad(n,a,b);
    for jj=1:length(quad_pts)
        term_quad(jj)=quad_wts(jj)*subs(integrand,r,quad_pts(jj));
    end
    result=vpa(sum(term_quad));
end

function [quad_pts,quad_wts]=giveGaussJacobiQuad(n,a,b)
    syms r
    jacobi_poly_n=jacobiP(n,a,b,r);
    quad_pts_0=vpasolve(jacobi_poly_n);
    quad_pts=(quad_pts_0+1)/2; 
    jacobi_poly_n_plus_1=jacobiP(n+1,a,b,r);
    for ii=1:length(quad_pts)
        quad_wts(ii)=-((2*n+a+b+2)/(n+a+b+1))*((gamma(n+a+1)*gamma(n+b+1))/(gamma(n+a+b+1)*factorial(n+1)))*(2^(a+b))/subs(diff(jacobi_poly_n,r)*jacobi_poly_n_plus_1,r,quad_pts_0(ii));
    end
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
    polys_to_use=sym("0");
    poly_coeff=sym("nu",[N N]);
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
    kk=1;
    poly_coeff_in_use=sym("0");
    for ii=2:N
        for jj=1:ii-1
            poly_coeff_in_use(kk)=poly_coeff(ii,jj);
            kk=kk+1;
        end
    end
    ortho_eqn=sym("0");
    kk=1;
    for ii=1:N
        for jj=1:N
            if(ii~=jj && ii<jj)           
                ortho_eqn(kk)=int(polys_to_use(ii)*polys_to_use(jj)*r^2,r,0,1)==0;
                kk=kk+1;
            end
        end
    end
    assume(poly_coeff,"real");
    poly_coeff_sol=solve(ortho_eqn);
    if(N>2)  
        poly_coeff_sol=cell2sym(struct2cell(poly_coeff_sol))';
    end
    polys_to_use=subs(polys_to_use,poly_coeff_in_use,poly_coeff_sol);
    SP_polys_to_use=polys_to_use;
    thiag_polys_to_use=SP_polys_to_use(2:end);
    rootP=vpasolve(thiag_polys_to_use(end));
    rootP=sort(abs(rootP(1:end/2)));
    rootP=[rootP' 1];
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
        quad_wts(ii)=int(Lgr_wts(ii)*3*r^2,r,0,1);
    end
end

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

%%

[Capprox_SP, ~]=giveCapproxNew(N);
[Capprox_thiag, ~]=giveCapproxThiag(N);
[Capprox_SP_high, ~]=giveCapproxNew(N_high);

Cs_sol_SP=zeros(1,length(taupoints))';
Cs_sol_SP_mix=zeros(1,length(taupoints))';

Cavg_sol_SP=zeros(1,length(taupoints));
Cavg_sol_SP_mix=zeros(1,length(taupoints));
Cavg_thiag_new=zeros(1,length(taupoints));
Cavg_thiag_GJ=zeros(1,length(taupoints));

SP_s = matlabFunction(subs(subs(Capprox_SP,r,1),alpha,aa));
SP_s_mix = matlabFunction(subs(subs(Capprox_SP_high,r,1),alpha_high,aa_high));

SP_s_high = matlabFunction(subs(subs(Capprox_SP_high,r,1),alpha_high,aa_high));
SP_s_low = matlabFunction(subs(subs(Capprox_SP_low,r,1),alpha_low,aa_low));

Cavg_thiag_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));
Cavg_thiag_GJ_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));

Cs_thiag_new=c_sol_thiag(:,end);
Cs_thiag_GJ=c_sol_thiag_GJ(:,end);

alpha_sol_SP_mix_orig=alpha_sol_SP_mix;

for tt=1:length(taupoints)
    alpha_for_SP=num2cell(alpha_sol_SP(tt,:));
    Cs_sol_SP(tt)=SP_s(alpha_for_SP{:});

    alpha_for_SP_mix=num2cell(alpha_sol_SP_mix(tt,:));
    Cs_sol_SP_mix(tt)=SP_s_mix(alpha_for_SP_mix{:});

    Cavg_sol_SP(tt)=alpha_sol_SP(tt,1);
    Cavg_sol_SP_mix(tt)=alpha_sol_SP_mix(tt,1);

    c_for_thiag=num2cell(c_sol_thiag(tt,:));
    Cavg_thiag_new(tt)=Cavg_thiag_fn(c_for_thiag{:});

    c_for_thiag_GJ=num2cell(c_sol_thiag_GJ(tt,:));
    Cavg_thiag_GJ(tt)=Cavg_thiag_GJ_fn(c_for_thiag_GJ{:});
end

alpha_sol_SP_mix_temp=alpha_sol_SP_mix_orig;


for tt=1:length(alpha_sol_total)
    alpha_for_mix_new=num2cell(alpha_sol_total(tt,:));
    Cs_sol_mix_new(tt)=Csurf_high_fun(alpha_for_mix_new{:});
    Cavg_sol_mix_new(tt)=alpha_sol_total(tt,1);
    flux_at_tspan_run(tt) = giveDeltafortau(tspan_run(tt), pStruct.Deltavstauprofile, pStruct);
    trig_switch_val(tt)=abs(-flux_at_tspan_run(tt)/ff_high(alpha_for_mix_new{:}));
end

% figure

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


% abs(-flux/ff_high(alpha_sol_SP_mix{:}))%% Plot the difference between Models

model_names=["\bfPDEPE\rm-\bf"+length(pStruct.rpoints)+"\rmnodes",...
    "\bf"+N+"\rmEqn"+"-\bfThiagOriginal\rm-\bf"+(N-1)+"\rmQdPts",...
    "\bf"+N+"\rmEqn"+"-\bfThiagModified\rm-\bf"+N_quad_thiag+"\rmQdPts",...
    "\bf"+N+"\rmEqn"+"-\bfNew\rm-\bf"+N_quad+"\rmQdPts",...
    "\bf("+N_low+"+"+N_high+")\rmEqn"+"-\bfMixNew\rm-\bf("+N_quad_low+"+"+N_quad_high+")\rmQdPts",
    ];


Cs_All= [C_PDEPE(:,end), Cs_thiag_new, Cs_thiag_GJ, Cs_sol_SP];
Cavg_All = [Cavg_PDEPE; Cavg_thiag_new; Cavg_thiag_GJ; Cavg_sol_SP];

clr_used=["blue", "black", "green", "magenta", "cyan"];

[~, cols]=size(Cs_All);
num_models=cols;
% 
% figure
% sgtitle("Comparison of N = "+N+" Models", 'FontSize',24,'FontWeight', 'bold');
% 
% taus_idx = [ceil((1/5) * length(taupoints)), 2*377, ...
%             2*544, ceil((4/5) * length(taupoints))];
% 
% profile_plt_loc = [15 16 21 22];
% lnWdt=1.5;
% for idx = 1:4
%     subplot(num_models, num_models + 2, profile_plt_loc(idx));
%     plot(rpoints, C_PDEPE(taus_idx(idx), :), "b",'LineWidth', lnWdt);
%     hold on;
%     C_thiag_new_vs_r = double(subs(subs(Capprox_thiag, c, c_sol_thiag(taus_idx(idx), :)), r, rpoints));
%     plot(rpoints, C_thiag_new_vs_r, "--k",'LineWidth', lnWdt);
%     hold on;
%     C_thiag_GJ_vs_r = double(subs(subs(Capprox_thiag, c, c_sol_thiag_GJ(taus_idx(idx), :)), r, rpoints));
%     plot(rpoints, C_thiag_GJ_vs_r, "-.g",'LineWidth', lnWdt);
%     hold on;
%     C_SP_vs_r = double(subs(subs(Capprox_SP, alpha, alpha_sol_SP(taus_idx(idx), :)), r, rpoints));
%     plot(rpoints, C_SP_vs_r, "--m",'LineWidth', lnWdt);
%     hold on;
%     C_mix_vs_r = double(subs(subs(Capprox_high, alpha_high, alpha_sol_total(taus_idx(idx), :)), r, rpoints));
%     plot(rpoints, C_mix_vs_r, "--c",'LineWidth', lnWdt);
% 
%     ylim([0 1.5]);
%     xlabel("r", FontSize = 16);
%     ylabel("C", FontSize = 16);
%     title("C vs r (\tau= " + taupoints(taus_idx(idx))+")");
% end
% 
% 
% % Calculate the number of unique pairs (N choose 2)
% num_pairs = nchoosek(num_models, 2);
% pair_idx = 1;
% err_plt_loc = pair_idx;
% 
% 
% % Iterate over unique pairs of models
% for ii = 1:num_models
%     for jj = ii+1:num_models
%         % Calculate the position for the subplot
%         row_shift = ceil((err_plt_loc + 2) / (num_models + 2));
% 
%         % Adjust the subplot position with the calculated shift
%         subplot(num_models, num_models + 2, err_plt_loc + 2 * row_shift);
% 
%         % Plot the surface concentration error on the left axis
%         yyaxis left
%         plot(taupoints, Cs_All(:, ii) - Cs_All(:, jj), 'LineWidth', 1,'Marker', 'none');
%         ylabel("\epsilon_{C_{s}}", 'FontSize', 12);
%         hold on
% 
%         % Plot the average concentration error on the right axis
%         yyaxis right
%         plot(taupoints, Cavg_All(ii, :) - Cavg_All(jj, :), 'LineWidth', 1 ,'Marker', 'none');
%         ylabel("\epsilon_{C_{avg}}", 'FontSize', 12);
% 
%         % Add error sum text for both errors
%         Cs_err_sqr_sum = sum((Cs_All(:, ii) - Cs_All(:, jj)).^2);
%         Cavg_err_sqr_sum = sum((Cavg_All(ii, :) - Cavg_All(jj, :)).^2);
%         Cs_err_text = "$\sum (\epsilon_{C_{\mathrm{s}}})^2: $" + sprintf('%.2e', Cs_err_sqr_sum);
%         text(0.2, 0.1, Cs_err_text, 'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'Latex', 'VerticalAlignment', 'top');
% 
%         Cavg_err_text="$\sum (\epsilon_{C_{\mathrm{avg}}})^2: $" + sprintf('%.2e', Cavg_err_sqr_sum);
%         text(0.2, 0.9, Cavg_err_text, 'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'Latex', 'VerticalAlignment', 'top');
% 
% 
%         % Add title and labels
%         xlabel("\tau", 'FontSize', 12);
%         title("\color{" + clr_used(ii) + "}[" + model_names(ii) + "]\color{black} - " + ...
%               "\color{" + clr_used(jj) + "}[" + model_names(jj) + "]", 'FontSize', 8);
% 
%         % Increment the error plot location
%         err_plt_loc = err_plt_loc + 1;
%     end
% end
% 
% subplot(num_models,num_models+2,[7 8 13 14])
% 
% yyaxis right
% % Plot Flux vs Time in the background
% flux_at_tau = zeros(1, length(taupoints));
% for ii = 1:length(taupoints)
%     flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
% end
% plot(taupoints, flux_at_tau, 'k:', 'LineWidth', 1.1, 'Color', [0 0 0 0.8]) % Dashed black line with transparency
% yticklabels([])
% ax = gca;
% ax.YColor = 'k';
% 
% yyaxis left
% plot(taupoints,C_PDEPE(:,end),"b",'LineWidth', lnWdt)
% hold on
% plot(taupoints,Cs_thiag_new,"--k",'LineWidth', lnWdt)
% hold on
% plot(taupoints,Cs_thiag_GJ,"-.g",'LineWidth', lnWdt)
% hold on
% plot(taupoints,Cs_sol_SP,"--m",'LineWidth', lnWdt)
% hold on
% % plot(taupoints,Cs_sol_SP_mix,"-.c",'LineWidth', lnWdt)
% plot(tspan_run,Cs_sol_mix_new,"-.c",'LineWidth', lnWdt)
% hold on
% xline(te_list,"r--")
% 
% 
% 
% for ii=1:length(taus_idx)
%     hold on
%     xline(taupoints(taus_idx(ii)),"--",LineWidth=0.5)
% end
% xlabel("\tau",FontSize=16)
% ylabel("C_s",FontSize=16)
% title("Surface Conc. vs time")
% legend_handle = legend(model_names, 'Location', 'northwest');
% legend_children = findobj(legend_handle, 'type', 'line');
% ax = gca;
% ax.YColor = 'k';
% 
% subplot(num_models,num_models+2,19)
% model_runtimes = [PDEPE_time, thiag_time, thiag_GJ_time, SP_time, Mix_New_time];
% 
% bar_handle = barh(model_runtimes, 'FaceColor', 'flat');
% yticklabels(model_names)
% 
% 
% bar_handle.CData = [
%     0, 0, 1;       % blue
%     0, 0, 0;       % black
%     0, 1, 0;       % green
%     1, 0, 1;       % magenta
%     0, 1, 1;       % cyan
% ];
% 
% set(gca,'XScale', 'log');
% title('Model Runtimes', 'FontSize', 10);
% 
% xlim([1e-3 1e2]); % Adjust x-axis limits for better visibility
% 
% % Add runtime values next to each bar
% for ii = 1:length(model_runtimes)
%     text(model_runtimes(ii), ii, sprintf('%.2f ms', model_runtimes(ii) * 1000), ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 8);
% end
% 
% set(gca,'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
% 
% 
% subplot(num_models,num_models+2,1)
% semilogy(CC,f_at_C)
% xlabel("C",FontSize=16)
% ylabel("f",FontSize=16)
% title("Diffusivity vs Conc.")
% 
% flux_at_tau=zeros(1,length(taupoints));
% for ii=1:length(taupoints)
%     flux_at_tau(ii)=giveDeltafortau(taupoints(ii),pStruct.Deltavstauprofile,pStruct);
% end
% 
% subplot(num_models,num_models+2,2)
% plot(taupoints,flux_at_tau)
% xlabel("\tau",FontSize=16)
% ylabel("\delta",FontSize=16)
% title("Flux vs Time")
% 
% subplot(num_models,num_models+2,20)
% plot(taupoints,Cavg_PDEPE,"b")
% hold on
% plot(taupoints,Cavg_thiag_new,"--k")
% hold on
% plot(taupoints,Cavg_thiag_GJ,"-.g")
% hold on
% plot(taupoints,Cavg_sol_SP,"--m")
% hold on
% % plot(taupoints,Cavg_sol_SP_mix,"-.c")
% plot(tspan_run,Cavg_sol_mix_new,"-.c")
% 
% 
% for ii=1:length(taus_idx)
%     hold on
%     xline(taupoints(taus_idx(ii)),"--",LineWidth=0.5)
% end
% xlabel("\tau",FontSize=16)
% ylabel("C_{avg}",FontSize=16)
% title("Average Conc. vs Time")
% 
% subplot(num_models,num_models+2,[11 12])
% yyaxis right
% plot(taupoints, flux_at_tau, 'k:', 'LineWidth', 2, 'Color', [0 0 0 0.3]) % Dashed black line with transparency
% ylabel("\delta", 'FontSize', 16, 'Color', [0.3 0.3 0.3]) % Grey color for text
% title("Modes vs Time"+" \color{" + clr_used(4) + "}[" + model_names(4) + "]")
% ax = gca;
% ax.YColor = 'k';
% 
% % Plot Modes vs Time in the foreground
% yyaxis left
% vibgyor_colors = [1 0 0; 1 0.5 0; 0.85 0.65 0; 0 1 0; 0 0 1; 0.29 0 0.51; 0.56 0 1; 0.5 0.5 0.5; 0 0.75 0.75; 0.75 0 0.75]; % Extended spectrum
% num_colors = size(vibgyor_colors, 1);
% for ii = 1:size(alpha_sol_SP, 2)
%     color_idx = num_colors - mod(ii-1, num_colors); % Reverse VIBGYOR order
%     plot(taupoints, alpha_sol_SP(:, ii), 'LineStyle', '-', 'LineWidth', 1, ...
%          'Color', vibgyor_colors(color_idx, :), 'Marker', 'none') % Ensure no markers
%     hold on
% end
% xlabel("\tau", 'FontSize', 16)
% ylabel("\alpha", 'FontSize', 16)
% legend(arrayfun(@(x) strcat("\alpha_", num2str(x)), 1:size(alpha_sol_SP, 2), 'UniformOutput', false), 'Location', 'best')
% ax = gca;
% ax.YColor = 'k';
% 
% 
% subplot(num_models,num_models+2,[17 18])
% flux_at_tau = zeros(1, length(taupoints));
% for ii = 1:length(taupoints)
%     flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
% end
% 
% yyaxis right
% plot(taupoints, flux_at_tau, 'k:', 'LineWidth', 2, 'Color', [0 0 0 0.3]) % Dashed black line with transparency
% ylabel("\delta", 'FontSize', 16, 'Color', [0.3 0.3 0.3]) % Grey color for text
% title("Ratio of Modes vs Time"+" \color{" + clr_used(4) + "}[" + model_names(4) + "]")
% ax = gca;
% ax.YColor = 'k';
% 
% % Plot Modes vs Time in the foreground
% yyaxis left
% num_colors = size(vibgyor_colors, 1);
% for ii = 2:size(alpha_sol_SP, 2)
%     color_idx = num_colors - mod(ii-1, num_colors); % Reverse VIBGYOR order
%     alpha_ratio=alpha_sol_SP(:, ii)./alpha_sol_SP(:,ii-1);
%     plot(taupoints,alpha_ratio , 'LineStyle', '-', 'LineWidth', 1, ...
%          'Color', vibgyor_colors(color_idx, :), 'Marker', 'none') % Ensure no markers
%     hold on
% end
% xlabel("\tau", 'FontSize', 16)
% ylabel("\alpha_{i}/\alpha_{i-1}", 'FontSize', 16)
% legend(arrayfun(@(x) strcat("\alpha_", num2str(x),"/\alpha_", num2str(x-1)), 2:size(alpha_sol_SP, 2), 'UniformOutput', false), 'Location', 'best')
% ax = gca;
% ax.YColor = 'k';
% 
% 
% subplot(num_models,num_models+2,[23 24])
% power_coeff_SP = zeros(length(taupoints),size(alpha_sol_SP,2)); 
% 
% Coeff_SP_poly = coeffs(Capprox_SP,r);
% Coeff_SP_poly = subs(Coeff_SP_poly,alpha,aa);
% Coeff_SP_fn=matlabFunction(Coeff_SP_poly);
% 
% for ii=1:length(taupoints)
%     alpha_temp=num2cell(alpha_sol_SP(ii,:));
%     power_coeff_SP(ii,:) = Coeff_SP_fn(alpha_temp{:});
% end
% 
% flux_at_tau = zeros(1, length(taupoints));
% for ii = 1:length(taupoints)
%     flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
% end
% 
% yyaxis right
% plot(taupoints, flux_at_tau, 'k:', 'LineWidth', 2, 'Color', [0 0 0 0.3]) % Dashed black line with transparency
% ylabel("\delta", 'FontSize', 16, 'Color', [0.3 0.3 0.3]) % Grey color for text
% title("Coeff of Powers vs Time"+" \color{" + clr_used(4) + "}[" + model_names(4) + "]")
% ax = gca;
% ax.YColor = 'k';
% 
% yyaxis left
% num_colors = size(vibgyor_colors, 1);
% for ii = 1:size(power_coeff_SP, 2)
%     color_idx = num_colors - mod(ii-1, num_colors); % Reverse VIBGYOR order
%     plot(taupoints,power_coeff_SP(:,ii), 'LineStyle', '-', 'LineWidth', 1, ...
%          'Color', vibgyor_colors(color_idx, :), 'Marker', 'none') % Ensure no markers
%     hold on
% end
% xlabel("\tau", 'FontSize', 16)
% ylabel("\beta", 'FontSize', 16)
% legend(arrayfun(@(x) strcat("\beta_", num2str(x)), 1:size(power_coeff_SP, 2), 'UniformOutput', false), 'Location', 'best')
% ax = gca;
% ax.YColor = 'k';

% figure;
label_fontsize = 8;
tick_fontsize = 12;
legend_fontsize = 12;
title_fontsize = 30; % New: font size for plot titles
fontname = 'Cambria Math';
interpreter = 'latex';
fig_width_cm = 16; % width in cm
fig_height_cm = 25; % height in cm (adjust as needed)
fig_width_in = fig_width_cm / 2.54; % convert to inches
fig_height_in = fig_height_cm / 2.54;

h_gap=0.13;
v_gap=0.07;
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

plot(taupoints,Cs_sol_SP,"--m",'LineWidth', 2)
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$C_s$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Surface~Conc.~vs~Time}$', 'Interpreter', interpreter); % Set title handle
set(th, 'FontSize', title_fontsize, 'FontName', fontname); % Explicitly set font size and font name

% Cavg_sol_JP=zeros(1,length(taupoints));
Cavg_sol_SP=zeros(1,length(taupoints));
Cavg_thiag_new=zeros(1,length(taupoints));
Cavg_thiag_GJ=zeros(1,length(taupoints));

%create a matlab function out of the symbolic experssion (for faster computation)
% JP_s = matlabFunction(subs(subs(Capprox_JP,r,1),alpha,aa));
SP_s = matlabFunction(subs(subs(Capprox_SP,r,1),alpha,aa));

Cavg_thiag_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));
Cavg_thiag_GJ_fn=matlabFunction(int(3*r^2*subs(Capprox_thiag,c,cc),r,0,1));

%for thiag its just the last c, cs=c(end)
Cs_thiag_new=c_sol_thiag(:,end);
Cs_thiag_GJ=c_sol_thiag_GJ(:,end);

%for SM the variables themselves are just Cavg and Cs
% Cavg_SM=c_sol_SM(:,1);
% Cs_SM=c_sol_SM(:,2);

for tt=1:length(taupoints)
    %surface conc.
    % alpha_for_JP=num2cell(alpha_sol_JP(tt,:));
    % Cs_sol_JP(tt)=JP_s(alpha_for_JP{:});
    alpha_for_SP=num2cell(alpha_sol_SP(tt,:));
    Cs_sol_SP(tt)=SP_s(alpha_for_SP{:});

    % average conc. equal to alpha_1
    % Cavg_sol_JP(tt)=alpha_sol_JP(tt,1);
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
    ['$\mathbf{$', num2str(N), '$}\;\mathrm{Eqn}-\mathbf{ThiagOriginal}-$', num2str(N-1), '$\;\mathrm{QdPts}$'], ...
    ['$\mathbf{$', num2str(N), '$}\;\mathrm{Eqn}-\mathbf{ThiagModified}-$', num2str(N_quad_thiag), '$\;\mathrm{QdPts}$'], ...
    ['$\mathbf{$', num2str(N), '$}\;\mathrm{Eqn}-\mathbf{New}-$', num2str(N_quad), '$\;\mathrm{QdPts}$'] ...
    };


Cs_sol_mix_interp=interp1(tspan_run,Cs_sol_mix_new,taupoints);
Cavg_sol_mix_interp=interp1(tspan_run,Cavg_sol_mix_new,taupoints);


Cs_All=[C_PDEPE(:,end), Cs_thiag_new, Cs_thiag_GJ, Cs_sol_SP, Cs_sol_mix_interp'];
Cavg_All = [Cavg_PDEPE; Cavg_thiag_new; Cavg_thiag_GJ; Cavg_sol_SP; Cavg_sol_mix_interp];

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
    "\bf"+(N)+"\rmN"+"-\bfThiagOriginal\rm-\bf"+(N-1)+"\rmQdPts",...
    "\bf"+(N)+"\rmN"+"-\bfThiagModified\rm-\bf"+N_quad_thiag+"\rmQdPts",...
    "\bf"+N+"\rmN"+"-\bfNewModal\rm-\bf"+N_quad+"\rmQdPts",...
    "\bf("+N_low+"+"+N_high+")\rmN"+"-\bfMixModal\rm-\bf("+N_quad_low+"+"+N_quad_high+")\rmQdPts"];


Cs_err_with_time=zeros(6,length(taupoints));
Cavg_err_with_time=zeros(6,length(taupoints));

% Cs_err_mix_with_time=Cs_All(:,1)-Cs_sol_mix_new(:;

for kk=[11,12]
    for ii = 1:num_models
        for jj = ii+1:num_models
 
            Cs_err_with_time(pair_idx,:)=Cs_All(:,ii)-Cs_All(:,jj);
            Cavg_err_with_time(pair_idx,:)=Cavg_All(ii,:)-Cavg_All(jj,:);
    
            pair_idx=pair_idx+1;
        end
    end
    subtightplot(6, 2, kk, [v_gap h_gap], [bt_mrg tp_mrg], [lf_mrg rt_mrg]);
    if(kk==11)
        plot(taupoints,Cs_err_with_time(1,:),'-k', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cs_err_with_time(2,:),'-g', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cs_err_with_time(3,:),'--m', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cs_err_with_time(4,:),'--c', 'LineWidth', 1 ,'Marker', 'none');
        xlabel('$\tau$', 'FontSize', label_fontsize-4, 'FontName', fontname, 'Interpreter', interpreter);
        ylabel('${C_{\mathrm{s,PDEPE}}}-{C_{\mathrm{s,Model}}}$', 'FontSize', label_fontsize-4, 'FontName', fontname, ...
            'Interpreter', interpreter);
        th = title("Error in Surface Conc."); % Set title handle
        set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
        set(gca, 'FontName', fontname, 'FontSize', tick_fontsize, 'TickLabelInterpreter', interpreter);
    else
        plot(taupoints,Cavg_err_with_time(1,:),'-k', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cavg_err_with_time(2,:),'-g', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cavg_err_with_time(3,:),'--m', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        plot(taupoints,Cavg_err_with_time(4,:),'--c', 'LineWidth', 1 ,'Marker', 'none');
        hold on;
        xlabel('$\tau$', 'FontSize', label_fontsize-4, 'FontName', fontname, 'Interpreter', interpreter);
        ylabel('${C_{\mathrm{avg,PDEPE}}}-{C_{\mathrm{avg,Model}}}$', 'FontSize', label_fontsize-4, 'FontName', fontname, ...
            'Interpreter', interpreter);
        th = title("Error in Average Conc."); % Set title handle
        set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
        set(gca, 'FontName', fontname, 'FontSize', tick_fontsize, 'TickLabelInterpreter', interpreter);
        % ax=gca;
        ax = gca;
        ax.YAxis.Exponent = -7;

    end
end

subtightplot(6,2,[7 8],[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
plot(tspan_run,trig_switch_val,"-r");
xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
% ylabel('|\frac{-\delta(\tau)}{f\(r=1\)}|', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)

ylabel('$\left| -\frac{\delta(\tau)}{f(C_s)} \right|$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
set(gca, 'FontSize', tick_fontsize, 'FontName', fontname);
hold on
xline(te_list,"--r")
th = title('$\mathrm{Surface~Conc.~Slope~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);



subtightplot(6,2,[3 4 5 6],[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])

yyaxis right
flux_at_tau = zeros(1, length(taupoints));
for ii = 1:length(taupoints)
    flux_at_tau(ii) = giveDeltafortau(taupoints(ii), pStruct.Deltavstauprofile, pStruct);
end
plot(taupoints, flux_at_tau, ':k', 'LineWidth', 1.1, 'Color', [0 0 0 0.8])
yticklabels([])
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
hold on
plot(tspan_run,Cs_sol_mix_new,"--c",'LineWidth', lnWdt)
hold on
xline(te_list,"r--")


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

subtightplot(6,2,10,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
model_runtimes = [PDEPE_time, thiag_time, thiag_GJ_time, SP_time, Mix_New_time];
bar_handle = bar(model_runtimes, 0.5,'FaceColor', 'flat');
set(gca, 'FontSize', tick_fontsize, 'FontName', fontname);

bar_handle.CData = [
    0, 0, 1;       % blue
    0, 0, 0;       % black
    0, 1, 0;       % green
    1, 0, 1;       % magenta
    0, 1, 1;       % cyan
];

set(gca,'YScale', 'log');
ylabel('$\mathrm{Runtime~[s]}$', 'FontSize', label_fontsize-6, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'TickLabelInterpreter', interpreter, 'FontSize', tick_fontsize, 'FontName', fontname);
th = title('$\mathrm{Model~Runtimes}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
ylim([1e-3 1e3]);
for ii = 1:length(model_runtimes)
    text(ii, model_runtimes(ii)*1.1, sprintf('$%.0f\\mathrm{ms}$', model_runtimes(ii) * 1000), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, ...
        'FontName', fontname, 'Interpreter', interpreter);
end
set(gca,'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'FontName', fontname, 'FontSize', tick_fontsize);
xticklabels([])


subtightplot(6,2,9,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
plot(taupoints,Cavg_PDEPE,"b")
hold on
plot(taupoints,Cavg_thiag_new,"--k")
hold on
plot(taupoints,Cavg_thiag_GJ,"-.g")
hold on
plot(taupoints,Cavg_sol_SP,"--m")
hold on
plot(tspan_run,Cavg_sol_mix_new,"--c")

xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$C_{\mathrm{avg}}$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Average~Conc.~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)

subtightplot(6,2,2,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
plot(taupoints,flux_at_tau)
xlabel('$\tau$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$\delta$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Flux~vs~Time}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)
ylim([-0.04 0.08])

subtightplot(6,2,1,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
semilogy(CC,f_at_C)
xlabel('$C$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
ylabel('$f$', 'FontSize', label_fontsize, 'FontName', fontname, 'Interpreter', interpreter)
th = title('$\mathrm{Diffusivity~vs~Conc.}$', 'Interpreter', interpreter);
set(th, 'FontSize', title_fontsize, 'FontName', fontname, 'Interpreter', interpreter);
set(gca, 'FontName', fontname, 'FontSize', tick_fontsize)

flux_at_tau=zeros(1,length(taupoints));
for ii=1:length(taupoints)
    flux_at_tau(ii)=giveDeltafortau(taupoints(ii),pStruct.Deltavstauprofile,pStruct);
end
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

    C_mix_vs_r = double(subs(subs(Capprox_high, alpha_high, alpha_sol_total(taus_idx(idx), :)), r, rpoints));
    plot(rpoints, C_mix_vs_r, "--c",'LineWidth', lnWdt);


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