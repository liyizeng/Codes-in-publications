% This is part of the orginal codes used in the following paper:
% http://www.molbiolcell.org/cgi/doi/10.1091/mbc.E22-10-0494
% On the role of myosin-induced actin depolymerization during cell migration
% If you have questions, feel free to contact Dr. Yizeng Li.

if Selection == 1
    if Parameter1 == 0
        N1 = 1;
    else
        N1 = 35;
    end
    if Parameter2 == 0
        N2 = 1;
    else
        N2 = 52;
    end
    
    if Parameter1 == 1
        JACTINF0 = logspace(-1,1,N1);
        Ym = JACTINF0;
    elseif Parameter1 == 2
        ETAST = logspace(-2,1,N1)*etast;
        Ym = ETAST*1d18/1d12;
    elseif Parameter1 == 3
        GAMMA = logspace(-0,3,N1)*gamma0;
        if JactinMode == 1
            GAMMA = logspace(-3.2,-1,N1);
        elseif JactinMode == 2
            GAMMA = logspace(-1,-4,N1);
        end
        Ym = GAMMA;
    elseif Parameter1 == 4
        KAD = logspace(-2,1,N1);
        Ym = KAD*1d9/1d6;
    elseif Parameter1 == 5
        DG = logspace(-2,1,N1)*dg;
        Ym = DG*1d9/1d6;
    elseif Parameter1 == 6
        ETA = logspace(-2,1,N1)*eta;
        Ym = ETA*1d18/1d12;
    elseif Parameter1 == 7
        DTC = logspace(-2,1,N1)*Dtc;
        Ym = DTC*1d12/1d18;
    elseif Parameter1 == 8
        KSIGMAN = logspace(1,3,N1);
        Ym = KSIGMAN;
    elseif Parameter1 == 9
        FEXTF = linspace(0,15,N1);
        Ym = FEXTF;
    end
    
    if Parameter2 == 1
        JACTINF0 = logspace(-2,1,N2)*Jactinf0;
        Xm = JACTINF0;
    elseif Parameter2 == 2
        ETAST = logspace(-0.7,2.3,N2)*etast;
        Xm = ETAST*1d18/1d12;
    elseif Parameter2 == 3
        if JactinMode == 1
            GAMMA = logspace(-3.2,-1,N2);
        elseif JactinMode == 2
            GAMMA = logspace(-3.5,-1.3,N2);
            GAMMA = logspace(-0,-4,N2);
        end
        Xm = GAMMA;
    elseif Parameter2 == 4
        KAD = logspace(-2.01,1.01,N2);
        Xm = KAD*1d9/1d6;
    elseif Parameter2 == 5
        DG = logspace(-2,1,N2)*dg;
        Xm = DG*1d9/1d6;
    elseif Parameter2 == 6
        ETA = logspace(-2,1,N2)*eta;
        Xm = ETA*1d18/1d12;
    elseif Parameter2 == 7
        DTC = logspace(-2,1,N2)*Dtc;
        Xm = DTC*1d12/1d18;
    elseif Parameter2 == 8
        KSIGMAN = logspace(-2,1,N2)*ksigman;
        Xm = KSIGMAN;
    elseif Parameter2 == 9
        FEXTF = linspace(0,15,N2);
        Xm = FEXTF;
    end
    
elseif Selection == 2
    if Parameter2 == 1
        Jactinf0 = JACTINF0(loop2);
    elseif Parameter2 == 2
        etast = ETAST(loop2);
    elseif Parameter2 == 3
        gamma0 = GAMMA(loop2);
    elseif Parameter2 == 4
        kad = KAD(loop2);
    elseif Parameter2 == 5
        dg = DG(loop2);
    elseif Parameter2 == 6
        eta = ETA(loop2);
    elseif Parameter2 == 7
        Dtc = DTC(loop2);
    elseif Parameter2 == 8
        ksigman = KSIGMAN(loop2);
    elseif Parameter2 == 9
        fextf = FEXTF(loop2);
    end
    
elseif Selection == 3
    if Parameter1 == 1
        Jactinf0 = JACTINF0(loop1);
    elseif Parameter1 == 2
        etast = ETAST(loop1);
    elseif Parameter1 == 3
        gamma0 = GAMMA(loop1);
    elseif Parameter1 == 4
        kad = KAD(loop1);
    elseif Parameter1 == 5
        dg = DG(loop1);
    elseif Parameter1 == 6
        eta = ETA(loop1)*ones(N,1);
    elseif Parameter1 == 7
        Dtc = DTC(loop1);
    elseif Parameter1 == 8
        ksigman = KSIGMAN(loop1);
    elseif Parameter1 == 9
        fextf = FEXTF(loop1);
    end
    
end