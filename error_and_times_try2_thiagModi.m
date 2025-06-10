clear
clc

dd=8;
digits=dd;

pStruct.FvsCprofile=char("nmc811AMIDRFitted5");
pStruct.Deltavstauprofile=char("testCrate");
pStruct.C0vsrprofile=char("const");
pStruct.DeltaIsZeroForSeconds=20;
pStruct.Compress=1; 
pStruct.R=2e-6;
pStruct.c0=51554;
pStruct.D0=2e-14;
pStruct.D0multiplier=1e-3;
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


PDEPE_timer=tic;
[C_PDEPE, t_PDEPE] = givePDEPESol(rpoints,taupoints,pStruct);
PDEPE_time=toc(PDEPE_timer);

N_list=2:8;
% N_list=[2 4];

% m_list=[2 3];
QdPts_list=[1:5 6:2:24];
% QdPts_list=1:16;

% QdPts_list=[1:5];

% QdPts_list=1:10;
%extra qdpts points

err_Cs=zeros(length(N_list),length(QdPts_list));
err_Cavg=zeros(length(N_list),length(QdPts_list));
err_Cint=zeros(length(N_list),length(QdPts_list));

run_time=zeros(length(N_list),length(QdPts_list));
temp_run_time=run_time;

NTry=5;
tic
for iTry=1:NTry
    for iN=1:length(N_list)
        % N_list
        for iQd=1:length(QdPts_list)
    
            N=N_list(iN);
            % disp("N="+N+" QdPts="+QdPts_list(iQd));
            % QdPts_list(iQd);
            % syms r t
            t=sym("t");
            r=sym("r");
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
            Capprox_SP=subs(Capprox,P,SP_polys_to_use);
            Capprox_SP;
            Coeff_SP_poly = coeffs(Capprox_SP,r);
            Coeff_SP_poly = subs(Coeff_SP_poly,alpha,aa);
            
            Nc=N-1;
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
            cc=sym("c",[1 Nc+1]);
            c= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'),t), cc, 'UniformOutput', false));
            Capprox_thiag=subs(Capprox,[P alpha],[Lgr_wts c]);
            
            ff = @(r,Capprox) firstOutputOnly(r,Capprox,pStruct);
     
            QdPts=QdPts_list(iQd);
            N_quad_SP=QdPts;
            N_quad_thiag_GJ=QdPts;

            a=0;
            b=2;
    
            dlt=sym("delta");
            delta= cell2sym(arrayfun(@(xEl) symfun(strcat(char(xEl), '(t)'), t), dlt, 'UniformOutput', false));
            
            W_SP=SP_polys_to_use;
            W_thiag=Lgr_wts;
            
            % Galerk_eqn_SP=sym("0");
            % Galerk_eqn_thiag=sym("0");
            Galerk_eqn_thiag_GJ=sym("0");
            
            for ii=1:N
                % Galerk_eqn_SP(ii)=int(W_SP(ii)*3*r^2*diff(Capprox_SP,t),r,0,1)==subs(W_SP(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_SP(ii),r)*3*ff(r,Capprox_SP)*diff(Capprox_SP,r),N_quad_SP,a,b);
                Galerk_eqn_thiag_GJ(ii)=int(W_thiag(ii)*3*r^2*diff(Capprox_thiag,t),r,0,1)==subs(W_thiag(ii)*3*r^2*(-delta),r,1)-giveQuadResult(diff(W_thiag(ii),r)*3*ff(r,Capprox_thiag)*diff(Capprox_thiag,r),N_quad_thiag_GJ,a,b);

            end
            
            digits=dd;
            % Galerk_eqn_SP=vpa(Galerk_eqn_SP);
            Galerk_eqn_thiag_GJ=vpa(Galerk_eqn_thiag_GJ);
    
            
            % [M_SP,F_SP]=massMatrixForm(Galerk_eqn_SP,alpha);
            % FF_SP=odeFunction(F_SP,alpha,delta);
            % M_SP=odeFunction(M_SP,alpha);
            % opt_SP=odeset("Mass",M_SP,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);
            
            [M_thiag_GJ,F_thiag_GJ]=massMatrixForm(Galerk_eqn_thiag_GJ,c);
            FF_thiag_GJ=odeFunction(F_thiag_GJ,c,delta);
            M_thiag_GJ=odeFunction(M_thiag_GJ,c);
            opt_thiag_GJ=odeset("Mass",M_thiag_GJ,"RelTol",pStruct.RelTol,"AbsTol",pStruct.AbsTol,"MaxStep",pStruct.MaxStep);



            % alpha0_SP=zeros(1,N);
            c0_thiag_GJ=zeros(1,N);

            
            for ii=1:N
                % alpha0_SP(ii)=double(int(pStruct.C0*SP_polys_to_use(ii)*r^2,r,0,1)/int(SP_polys_to_use(ii)*SP_polys_to_use(ii)*r^2,r,0,1));
                c0_thiag_GJ(ii)=pStruct.C0;

            end
                    
            flux = @(t) giveDeltafortau(t,pStruct.Deltavstauprofile,pStruct);
            
            fluxVal = @(t) giveDeltaVal(t,pStruct);  
            fluxDiff1 = @(t) giveDeltaDiff1(t,pStruct);
            fluxDiff2 = @(t) giveDeltaDiff2(t,pStruct);
            fluxDiff3 = @(t) giveDeltaDiff3(t,pStruct);
            
            % SP_timer=tic;
            % [~,alpha_sol_SP]=ode15s(@(t,Alpha) FF_SP(t,Alpha,flux(t)),pStruct.taupoints,alpha0_SP,opt_SP);
            % SP_time=toc(SP_timer);
            

            thiag_GJ_timer=tic;
            [~,c_sol_thiag_GJ]=ode15s(@(t,CC) FF_thiag_GJ(t,CC,flux(t)),pStruct.taupoints,c0_thiag_GJ,opt_thiag_GJ);
            thiag_GJ_time=toc(thiag_GJ_timer);

            % thiag_time=0;
            % thiag_GJ_time=0;
            % 
            % temp_run_time(iN,iQd,:)=[thiag_time thiag_GJ_time SP_time];
            temp_run_time(iN,iQd)=thiag_GJ_time;
    
            Cs_thiag_GJ=c_sol_thiag_GJ(:,end);

            % Cs_sol_SP=zeros(1,length(taupoints))';
            % Cint_sol_SP=zeros(length(taupoints),length(rpoints));
            % 
            % Cavg_sol_SP=zeros(1,length(taupoints));
    
            % SP_s = matlabFunction(subs(subs(Capprox_SP,r,1),alpha,aa));
            % SP_int = matlabFunction(subs(Capprox_SP,alpha,aa));
    
            for tt=1:length(taupoints)
                % alpha_for_SP=num2cell(alpha_sol_SP(tt,:));
                % Cs_sol_SP(tt)=SP_s(alpha_for_SP{:});
                % Cint_sol_SP(tt,:)=SP_int(alpha_for_SP{:},rpoints);
                % 
                % Cavg_sol_SP(tt)=alpha_sol_SP(tt,1);

                c_for_thiag_GJ=num2cell(c_sol_thiag_GJ(tt,1:end-1));
                % Cavg_thiag_GJ(tt)=Cavg_thiag_GJ_fn(c_for_thiag_GJ{:});
            end
    
            % Cavg_PDEPE=zeros(1,length(taupoints));
            % for ii=1:length(taupoints)
            %     Cr2=C_PDEPE(ii,:).*rpoints.^2;
            %     Cavg_PDEPE(ii)=(3)*trapz(rpoints,Cr2);
            % end
    
            % CC=linspace(0,1,100);
            % f_at_C=zeros(1,length(CC));
            % for ii=1:length(CC)
            %     f_at_C(ii)=giveFforC(CC(ii),pStruct.FvsCprofile,pStruct);
            % end
    
            Cs_err_sqr_sum = sum((C_PDEPE(:,end) - Cs_thiag_GJ(:)).^2);
            err_Cs(iN,iQd)=Cs_err_sqr_sum;
    
            % Cint_err_sqr_sum =sum(sum((C_PDEPE-Cint_sol_SP).^2,2));
            % err_Cint(iN,iQd)=Cint_err_sqr_sum;
        end
        if(iN==length(N_list))
                s = struct("run_time",temp_run_time,"err_Cs",err_Cs,"err_Cint",err_Cint);
                save(sprintf("err_times_data_"+"ThiagModi_"+pStruct.FvsCprofile+"_"+pStruct.Deltavstauprofile+"_Try"+iTry+".mat"),"-fromstruct",s)
        end
        % run_time=temp_run_time;
    end
end
toc

%function definitions
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


%% 
clear
clc

label_fontsize = 16;
tick_fontsize = 14;
legend_fontsize = 9;
fontname = 'Cambria Math';
interpreter = 'latex';
fig_width_cm = 16; % width in cm
fig_height_cm = 10; % height in cm (adjust as needed)
fig_width_in = fig_width_cm / 2.54; % convert to inches
fig_height_in = fig_height_cm / 2.54;


h_gap=0.15;
v_gap=0.11;
lf_mrg=0.15;
rt_mrg=0.01;
bt_mrg=0.15;
tp_mrg=0.1;

QdPts_list=[1:5 6:2:24];
% QdPts_list=1:16;
N_list=2:8;

%runtime plots with error bars
% rt_all_try=zeros(size(run_time,1),size(QdPts_list,2),5);
rt_all_try=[];
NTry=5;

for ii=1:NTry
    file_name="err_times_data"+"_NewModal_"+"nmc811AMIDRFitted5"+"_testCrate_Try"+ii;
    data1=load(file_name,"run_time");
    rt_all_try(:,:,ii)=data1.run_time;
end
data1=load(file_name,"err_Cs");
err_Cs=data1.err_Cs;

rt_mean=mean(rt_all_try,3);
rt_std_dev= std(rt_all_try,0,3);

f1 = figure();
set(f1, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
set(f1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 fig_width_in fig_height_in]);
set(f1, 'PaperSize', [fig_width_in fig_height_in]);

% Define a colormap from red to purple (using HSV, but only the red to purple segment)
nCurves = length(N_list);
% HSV: red (0), orange, yellow, green, blue, purple (0.8)
cmap = hsv2rgb([linspace(0,0.8,nCurves)' ones(nCurves,1) ones(nCurves,1)]);

subtightplot(1,2,1,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])

for ii=1:nCurves
    nm = "N = " + N_list(ii);

    % rt_mean_new=rt_mean(ii,ii:end);
    rt_mean_new=rt_mean(ii,:);
    % rt_std_dev_new=rt_std_dev(ii,ii:end);
    rt_std_dev_new=rt_std_dev(ii,:);
    % QdPts_list_new=QdPts_list(ii:end);
    QdPts_list_new=QdPts_list;
    errorbar(QdPts_list_new, rt_mean_new, rt_std_dev_new, 'DisplayName', nm, ...
        'LineWidth', 1, 'Color', cmap(ii,:));
    hold on
end
legend('Interpreter', interpreter)
title("Model Runtime")

ylabel("Runtime [s]",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
xlabel("Qd Pts",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
set(gca,'TickLabelInterpreter', interpreter)
set(gca, 'FontSize', tick_fontsize, 'FontName', fontname)

% --- Inset plot: zoomed in for x = 0 to 5 ---
axesPos = get(gca, 'Position'); % Get position of main axes
insetWidth = 0.35 * axesPos(3); % relative to main axes
insetHeight = 0.45 * axesPos(4);
insetX = axesPos(1) + 0.55 * axesPos(3);
insetY = axesPos(2) + 0.45 * axesPos(4);

insetAx = axes('Position', [insetX insetY insetWidth insetHeight]);
box on
hold on
for ii=1:nCurves
    errorbar(QdPts_list, rt_mean(ii,:), rt_std_dev(ii,:), ...
        'LineWidth', 1, 'Color', cmap(ii,:));
end
xlim([min(QdPts_list) 5])
set(insetAx, 'FontSize', tick_fontsize)
set(insetAx, 'FontName', fontname)
set(insetAx,'TickLabelInterpreter', interpreter)

% Move y-axis tick labels to the right
set(insetAx, 'YAxisLocation', 'right')

% Optionally, remove x/y labels for inset
% set(insetAx, 'XTickLabel', [])
% set(insetAx, 'YTickLabel', [])
hold off
conv_idx_list=[];


rt_all_try_thiagmodi=[];
NTry=5;

for ii=1:NTry
    file_name_thiagmodi="err_times_data"+"_ThiagModi_"+"nmc811AMIDRFitted5"+"_testCrate_Try"+ii;
    data1_thiagmodi=load(file_name_thiagmodi,"run_time");
    rt_all_try_thiagmodi(:,:,ii)=data1_thiagmodi.run_time;
end

data1_thiagmodi=load(file_name_thiagmodi,"err_Cs");
err_Cs_thiagmodi=data1_thiagmodi.err_Cs;

rt_mean_thiagmodi=mean(rt_all_try_thiagmodi,3);
rt_std_dev_thiagmodi= std(rt_all_try_thiagmodi,0,3);

figure
% subtightplot(1,2,2,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
for ii=1:nCurves
    nm = "N = " + N_list(ii);

    % rt_mean_new=rt_mean(ii,ii:end);
    rt_mean_new_thiagmodi=rt_mean_thiagmodi(ii,:);
    rt_mean_new=rt_mean(ii,:);

    % rt_std_dev_new=rt_std_dev(ii,ii:end);
    rt_std_dev_new_thiagmodi=rt_std_dev_thiagmodi(ii,:);
    % QdPts_list_new=QdPts_list(ii:end);
    QdPts_list_new_thiagmodi=QdPts_list;
    % errorbar(QdPts_list_new_thiagmodi, rt_mean_new_thiagmodi, rt_std_dev_new_thiagmodi, 'DisplayName', nm, ...
    %     'LineWidth', 1, 'Color', cmap(ii,:));
    plot(QdPts_list_new_thiagmodi,rt_mean_new-rt_mean_new_thiagmodi,'-o','MarkerSize',3,'DisplayName', nm, ...
   'LineWidth', 1, 'Color', cmap(ii,:))
    hold on
end
legend('Interpreter', interpreter)
title("Diff in Model Runtime (NewModal - ThiagModi)")

ylabel("Runtime [s]",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
xlabel("Qd Pts",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
set(gca,'TickLabelInterpreter', interpreter)
set(gca, 'FontSize', tick_fontsize, 'FontName', fontname)






% subtightplot(1,2,2,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
% for ii=1:nCurves
%     nm = "N = " + N_list(ii);
%     yvals = err_Cs(ii,:);
%     % yvals_new=yvals(ii:end);
%     yvals_new=yvals;
%     % QdPts_list_new=QdPts_list(ii:end);
%     QdPts_list_new=QdPts_list;
% 
%     semilogy(QdPts_list_new, yvals_new, 'DisplayName', nm, ...
%         'LineWidth', 1.5, 'Color', cmap(ii,:));
%     hold on
% 
%     % --- Find convergence point ---
%     % Define convergence as: relative change < 1% for 3 consecutive points
%     rel_change = abs(diff(yvals_new)) ./ abs(yvals_new(1:end-1));
%     % win = 3; % window of consecutive points
%     conv_idx = [];
%     for jj = 1:length(rel_change)
%         if rel_change(jj) < 0.01
%             conv_idx = jj+1; % mark the first point where convergence starts
%             break
%         end
%     end
%     if ~isempty(conv_idx)
%         % Mark with a star, but don't show in legend
%         plot(QdPts_list_new(conv_idx), yvals_new(conv_idx), 'p', ...
%             'MarkerSize', 10, 'MarkerEdgeColor', 'k', ...
%             'MarkerFaceColor', cmap(ii,:), 'HandleVisibility', 'off');
%     end
%     conv_idx_list=[conv_idx_list conv_idx];
% end
% legend('Interpreter', interpreter)
% title("Error in Surface Conc.")
% 
% ylabel("$\sum (\epsilon_{C_{\mathrm{s}}})^2$",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
% xlabel("Qd Pts",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
% set(gca,'TickLabelInterpreter', interpreter)
% set(gca, 'FontSize', tick_fontsize, 'FontName', fontname)

% subtightplot(1,3,3,[v_gap h_gap],[bt_mrg tp_mrg],[lf_mrg rt_mrg])
% 
% for ii=1:length(N_list)
%     nm="N = "+N_list(ii);
%     semilogy(QdPts_list,err_Cint(ii,:),'DisplayName', nm,'LineWidth',1.5)
%     hold on
% end
% legend
% % title("Error in Internal Conc. vs Number of Quadrature Points")
% title("Internal Conc.")
% 
% ylabel("$\sum (\epsilon_{C_{\mathrm{internal}}})^2$",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)
% xlabel("QdPts",'FontSize', label_fontsize, 'Interpreter', interpreter, 'FontName', fontname)

