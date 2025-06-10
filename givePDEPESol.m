function [C, timeElapsed] = givePDEPESol(r_mesh,tau_mesh,paramStruct,varargin)
    tic;
    FvsCprofile = paramStruct.FvsCprofile;
    Deltavstauprofile = paramStruct.Deltavstauprofile;
    C0vsrprofile = paramStruct.C0vsrprofile;

    function [c,f,s] = diffpde(x,t,u,dudx)
        [F,~] = giveFforC(u,FvsCprofile,paramStruct);
        c = 1;
        f = F*dudx;
        s = 0;
    end
    
    function u0 = diffic(x,~)
        u0 = giveC0forr(x,C0vsrprofile,paramStruct);
    end
    
    function [pl,ql,pr,qr] = diffbc(xl,ul,xr,ur,t)
        delta=giveDeltafortau(t,Deltavstauprofile,paramStruct);
        pl = 0;
        ql = 1;
        pr = delta;
        qr = 1;
    end

    m = 2; 
    RelTol=paramStruct.RelTol;
    AbsTol=paramStruct.AbsTol;

    MaxStep=paramStruct.MaxStep;

    opts=odeset("RelTol",RelTol,'AbsTol',AbsTol,'MaxStep',MaxStep);
    C = pdepe(m,@diffpde,@diffic,@diffbc,r_mesh,tau_mesh,opts);
    timeElapsed=toc;
end


