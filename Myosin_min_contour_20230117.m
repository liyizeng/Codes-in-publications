% This is part of the orginal codes used in the following paper:
% http://www.molbiolcell.org/cgi/doi/10.1091/mbc.E22-10-0494
% On the role of myosin-induced actin depolymerization during cell migration
% If you have questions, feel free to contact Dr. Yizeng Li.

% Two-phase: cytosol and actomyosin network
% Steady-state
% Unknowns:
% vn, theta_n, theta_c, v0
% The velocities are written in a fixed frame but the equations are written
% in the frame of the cell.

% 1: Jactinf
% 2: nust
% 4: kad
% 5: dg
% 6: eta
% 7: Dtc
% 8: ksigman
% 9: fextf

clear
clc

FadMode = 1;
% FadMode = 1: Fad = kad*v0
% FadMode = 2: Fad = kad*thetan(1)*v0

JactinMode = 2;
% JactinMode = 1: Jactinf = Jactinf0*thetac(N)
% JactinMode = 2: Jactinf = Jactinf0*thetac(N)/(thetacc + thetac(N))


RankCheck = 0;

N = 201;

Parameter1 = 1;
Parameter2 = 4;

%% Parameters are in units: nm, s, Pa & mM

L = 50.d3;              % (nm) cell length
p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell

Thetac = 0.1;           % (mM) reference value of G-actin
Thetan = 0.2;           % (mM) reference value of F-actin
Theta  = 0.3;           % (mM) reference total G- and F-actin, \int_0^L (thetan + thetac)dx/L

Mc = 5d-3;              % (mM) reference value of free myosin
Mn = 5d-3;              % (mM) reference value of

if JactinMode == 1
    Jactinf0 = 50;          % (nm/s) coefficient of actin polymerization
elseif JactinMode == 2
    Jactinf0 = 6;           % (nm mM/s) coefficient of actin polymerization
    thetacc  = 0.2d-3;         % (mM) Critical value for actin polymerization
end


ksigman = 100*1d2;          % (Pa /mM) Coefficient of passive actin network pressure
etast = 100*1d-4;           % (Pa s/nm^2/mM)
eta   = 1d-8;           % (Pa s/nm^2/mM)
dg    = 1d-6;           % (Pa s/nm) coefficient of hydraulic resistance
if FadMode == 1
    kad = 100*3d-1;           % (Pa s/nm) adhesive force, used in Lingxing's model
elseif FadMode == 2
    kad = 1d0;          % (Pa s/mM/nm) adhesive force, used in Lingxing's model
end
Dtc = 1.d7;             % (nm^2/s) diffusion constant for theta_c

%%
Iter = 21;

dx = L/(N-1);
x = linspace(0,L,N);

% number of variables of the model
n_var = 4; % v_n, theta_n, theta_c, v0
var_legh = [N,N,N,1]; % length of eavh variable
var_strt = zeros(1,n_var); % start position of each variable
var_strt(1) = 1;
for in = 2:n_var
    var_strt(in) = var_strt(in-1) + var_legh(in-1);
end
N_var = var_strt(n_var) + var_legh(n_var) - 1; % total size of the matrix

s_vn = var_strt(1); l_vn = var_legh(1);
s_tn = var_strt(2); l_tn = var_legh(2);
s_tc = var_strt(3); l_tc = var_legh(3);
s_v0 = var_strt(4); l_v0 = var_legh(4);

Selection = 1;
Sweep_20201003;

Ng = 51;
Gamma = logspace(-1,-4,Ng);          % (1/s)
V0 = zeros(Ng,1);
MAXGAMMA = zeros(N1,N2);

tic
Solution = 0;
for loop2 = 1:N2
    loop2
    toc
    Selection = 2;
    Sweep_20201003;
    for loop1 = 1:N1
        loop1;
        Selection = 3;
        Sweep_20201003;
        
        dgamma = 0.9;
        dv0 = inf;
        gamma0 = 1d-1; % (1/s)
        ig = 0;
        while dv0 > 1d-4
            ig = ig + 1;
            
            gamma = gamma0*logspace(-0,0,N)';
            nust = etast*logspace(-0,0,N)';
            
            % initial guess
            if (loop1 == 1 && loop2 == 1) || Solution == 0
                thetan = linspace(Thetan*1,Thetan*1,N)';
                thetac = linspace(Thetac*1,Thetac*1,N)';
                v0 = (1/2)*(etast*Theta*L)*(Jactinf0*gamma0*L)./((etast*Theta*L)*Jactinf0 ...
                    + (kad+dg)*(Jactinf0+gamma0*L));
                vn = linspace(v0,v0-gamma0*L,N)';
            elseif loop1 == 1 && loop2 > 1 && ig == 1
                vn = temp_vn_loop1;
                thetan = temp_thetan_loop1;
                thetac = temp_thetac_loop1;
                v0 = temp_v0_loop1;
            elseif ig == 1
                vn = temp_vn_ig;
                thetan = temp_thetan_ig;
                thetac = temp_thetac_ig;
                v0 = temp_v0_ig;
            end
            X = [vn; thetan; thetac; v0];
            
            iter = 0;
            ITER = true;
            while ITER
                iter = iter + 1;
                
                DF = zeros(N_var,N_var);
                Fn = zeros(N_var,1);
                
                if JactinMode == 1
                    Jactinf = Jactinf0*thetac(N);
                    DJactinfDthcN = Jactinf0;
                elseif JactinMode == 2
                    Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
                    DJactinfDthcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
                end
                
                %% Equations for vn and the Derivatives
                Fn(s_vn) = - ksigman*(thetan(2)-thetan(1)) + dx*eta*thetan(1)*(v0-vn(1)) ...
                    - dx*nust(1)*thetan(1)*vn(1);
                DF(s_vn,s_v0) = dx*eta*thetan(1);
                DF(s_vn,s_vn) = -dx*eta*thetan(1) - dx*nust(1)*thetan(1);
                DF(s_vn,s_tn) = ksigman + dx*eta*(v0-vn(1)) - dx*nust(1)*vn(1);
                DF(s_vn,s_tn+1) = -ksigman;
                
                Fn(s_vn+1:s_vn+N-2) = -ksigman*(thetan(3:N)-thetan(1:N-2)) + 2*dx*eta*thetan(2:N-1).*(v0-vn(2:N-1)) ...
                    - 2*dx*nust(2:N-1).*thetan(2:N-1).*vn(2:N-1);
                for i = 2:N-1
                    DF(s_vn+i-1,s_v0) = 2*dx*eta*thetan(i);
                    DF(s_vn+i-1,s_vn+i-1) = -2*dx*thetan(i)*(eta+nust(i));
                    DF(s_vn+i-1,s_tn+i-2) = ksigman;
                    DF(s_vn+i-1,s_tn+i-1) = 2*dx*eta*(v0-vn(i)) - 2*dx*nust(i)*vn(i);
                    DF(s_vn+i-1,s_tn+i)   = -ksigman;
                end
                
                Fn(s_vn+N-1) = thetan(N)*(vn(N) - v0) + Jactinf;
                DF(s_vn+N-1,s_vn+N-1) = thetan(N);
                DF(s_vn+N-1,s_tn+N-1) = vn(N) - v0;
                DF(s_vn+N-1,s_tc+N-1) = DJactinfDthcN;
                DF(s_vn+N-1,s_v0) = -thetan(N);
                
                %% Equations for thetan and the Derivatives
                Fn(s_tn) = thetan(1)*(vn(1)-v0);
                DF(s_tn,s_vn) = thetan(1);
                DF(s_tn,s_tn) = vn(1)-v0;
                DF(s_tn,s_v0) = -thetan(1);
                
                Fn(s_tn+1:s_tn+N-2) = (thetan(3:N).*vn(3:N)-thetan(1:N-2).*vn(1:N-2))...
                    - v0*(thetan(3:N)-thetan(1:N-2)) + 2*dx*gamma(2:N-1).*thetan(2:N-1);
                for i = 2:N-1
                    DF(s_tn+i-1,[s_vn+i-2,s_vn+i]) = [-thetan(i-1), thetan(i+1)];
                    DF(s_tn+i-1,s_tn+i-2) = -vn(i-1) + v0;
                    DF(s_tn+i-1,s_tn+i-1) = 2*dx*gamma(i);
                    DF(s_tn+i-1,s_tn+i)   =  vn(i+1) - v0;
                    DF(s_tn+i-1,s_v0) = -(thetan(i+1) - thetan(i-1));
                end
                
                Fn(s_tn+N-1) = (thetan(N)*vn(N)-thetan(N-1)*vn(N-1)) ...
                    - v0*(thetan(N)-thetan(N-1)) + dx*gamma(N)*thetan(N);
                DF(s_tn+N-1,[s_vn+N-2,s_vn+N-1]) = [-thetan(N-1), thetan(N)];
                DF(s_tn+N-1,s_tn+N-2) = -vn(N-1) + v0;
                DF(s_tn+N-1,s_tn+N-1) =  vn(N)   - v0 + dx*gamma(N);
                DF(s_tn+N-1,s_v0) = -(thetan(N)-thetan(N-1));
                
                %% Equations for thetac and the Derivatives
                Fn(s_tc) = thetac(2) - thetac(1);
                DF(s_tc,s_tc)   = -1;
                DF(s_tc,s_tc+1) =  1;
                
                Fn(s_tc+1:s_tc+N-2) = thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N)...
                    + dx^2/Dtc*gamma(2:N-1).*thetan(2:N-1);
                for i = 2:N-1
                    DF(s_tc+i-1,s_tn+i-1) = dx^2/Dtc*gamma(i);
                    DF(s_tc+i-1,s_tc+i-2) = 1;
                    DF(s_tc+i-1,s_tc+i-1) = -2;
                    DF(s_tc+i-1,s_tc+i)   = 1;
                end
                
                Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) + sum(thetan(2:N)+thetac(2:N)))...
                    - L*(Thetan + Thetac)/dx;
                DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
                DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
                DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
                DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;
                
                %% Equation for v0 and the Derivatives
                if FadMode == 1
                    Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad)*v0 ...
                        + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                    DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                    DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                    DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                    DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
                    DF(s_v0,s_v0) = dg + kad;
                elseif FadMode == 2
                    Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad*thetac(1))*v0 ...
                        + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                    DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                    DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                    DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                    DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
                    DF(s_v0,s_tn+1) = kad*v0;
                    DF(s_v0,s_v0) = dg + kad*thetac(1);
                end
                
                temp_Fn = Fn;
                
                %% Solve for the matrix
                DF = sparse(DF);
                X = X - DF\Fn;
                
                if sum(isnan(X)) >= 1 || iter == Iter || norm(temp_Fn) < norm(Fn)
                    Solution = 0;
                    vn     = NaN(N,1);
                    thetan = NaN(N,1);
                    thetac = NaN(N,1);
                    mn     = NaN(N,1);
                    mc     = NaN(N,1);
                    v0     = NaN;
                    break
                else
                    Solution = 1;
                end
                
                vn = X(s_vn:s_vn+N-1);
                thetan = X(s_tn:s_tn+N-1);
                thetac = X(s_tc:s_tc+N-1);
                v0 = X(s_v0);
                
                if iter > 1
                    error = abs((X-temp_X)./X);
                    error = sum(error)/(N_var);
                    if error < 1d-5 || iter == Iter
                        ITER = false;
                    end
                end
                temp_X = X;
            end
            
            if loop1 == 1 && ig == 1
                temp_vn_loop1 = vn;
                temp_thetan_loop1 = thetan;
                temp_thetac_loop1 = thetac;
                temp_v0_loop1 = v0;
            end
            if ig == 1
                temp_vn_ig = vn;
                temp_thetan_ig = thetan;
                temp_thetac_ig = thetac;
                temp_v0_ig = v0;
            end
            
            if ig == 1
                v1 = v0;
                g1 = gamma0;
                gamma0 = dgamma*gamma0;
            elseif ig == 2
                v2 = v0;
                g2 = gamma0;
                gamma0 = dgamma*gamma0;
            elseif ig == 3
                v3 = v0;
                g3 = gamma0;
                
                if v3>v2 && v2>v1
                    g1 = g2;
                    v1 = v2;
                    g2 = g3;
                    v2 = v3;
                    
                    gamma0 = gamma0*dgamma;
                else
                    glow  = g3;
                    ghigh = g1;
                    dv0 = 0;
                end
            end
            
            if ig > 3
                v_temp = v0;
                g_temp = gamma;
                v3 = v0;
                g3 = gamma0;
                if v0>v2 && v2>v1
                    g1 = g2;
                    v1 = v2;
                    g2 = g3;
                    v2 = v3;
                    gamma0 = gamma0*dgamma;
                    g3 = g_temp;
                    v3 = v_temp;
                else
                    glow  = g3;
                    ghigh = g1;
                    dv0 = 0;
                end
            end
                
                
                
                
%             if ig == 1     
%                 v0_temp1 = v0;
%                 gamma0_temp1 = gamma0;
%                 gamma0_temp = gamma0*(1-dgamma);
%             elseif ig == 2
%                 v0_temp2 = v0;
%                 gamma0_temp2 = gamma0;
%                 dv0 = (v0_temp2 - v0_temp1)/max(v0_temp1,v0_temp2);
%                 gamma0_temp = gamma0*(1+dgamma*(v0_temp2 - v0_temp1)/(gamma0_temp2-gamma0_temp1));
%             else
%                 v0_temp2 = v0;
%                 gamma0_temp2 = gamma0;
%                 dv0 = (v0_temp2 - v0_temp1)/max(v0_temp1,v0_temp2);
%                 gamma0_temp = gamma0*(1+dgamma*(v0_temp2 - v0_temp1)/(gamma0_temp2-gamma0_temp1));
%             end
%             
%             if ig > 1
%                 if gamma0_temp < 0
%                     ik = 1;
%                     while gamma0_temp < 0
%                         ik = ik + 1;
%                         gamma0_temp = gamma0*(1+dgamma/ik*(v0_temp2 - v0_temp1)/(gamma0_temp2-gamma0_temp1));
%                     end
%                 end
%                 v0_temp1 = v0;
%                 gamma0_temp1 = gamma0;
%                 
%             end
%             gamma0 = gamma0_temp;
%             
%             if Solution == 0
%                 gamma0_temp1 = NaN;
%             end

        end
       
        Gamma = linspace(glow,ghigh,Ng);
        for ig = 1:Ng
            gamma0 = Gamma(ig);
            
            gamma = gamma0*logspace(-0,0,N)';
            nust = etast*logspace(-0,0,N)';
            
            % initial guess
            if (loop1 == 1 && loop2 == 1 && ig == 1) || Solution == 0
                thetan = linspace(Thetan*1,Thetan*1,N)';
                thetac = linspace(Thetac*1,Thetac*1,N)';
                v0 = (1/2)*(etast*Theta*L)*(Jactinf0*gamma0*L)./((etast*Theta*L)*Jactinf0 ...
                    + (kad+dg)*(Jactinf0+gamma0*L));
                vn = linspace(v0,v0-gamma0*L,N)';
            elseif loop1 == 1 && loop2 > 1 && ig == 1
                vn = temp_vn_loop1;
                thetan = temp_thetan_loop1;
                thetac = temp_thetac_loop1;
                v0 = temp_v0_loop1;
            elseif ig == 1
                vn = temp_vn_ig;
                thetan = temp_thetan_ig;
                thetac = temp_thetac_ig;
                v0 = temp_v0_ig;
            end
            X = [vn; thetan; thetac; v0];
            
            iter = 0;
            ITER = true;
            while ITER
                iter = iter + 1;
                
                DF = zeros(N_var,N_var);
                Fn = zeros(N_var,1);
                
                if JactinMode == 1
                    Jactinf = Jactinf0*thetac(N);
                    DJactinfDthcN = Jactinf0;
                elseif JactinMode == 2
                    Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
                    DJactinfDthcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
                end
                
                %% Equations for vn and the Derivatives
                Fn(s_vn) = - ksigman*(thetan(2)-thetan(1)) + dx*eta*thetan(1)*(v0-vn(1)) ...
                    - dx*nust(1)*thetan(1)*vn(1);
                DF(s_vn,s_v0) = dx*eta*thetan(1);
                DF(s_vn,s_vn) = -dx*eta*thetan(1) - dx*nust(1)*thetan(1);
                DF(s_vn,s_tn) = ksigman + dx*eta*(v0-vn(1)) - dx*nust(1)*vn(1);
                DF(s_vn,s_tn+1) = -ksigman;
                
                Fn(s_vn+1:s_vn+N-2) = -ksigman*(thetan(3:N)-thetan(1:N-2)) + 2*dx*eta*thetan(2:N-1).*(v0-vn(2:N-1)) ...
                    - 2*dx*nust(2:N-1).*thetan(2:N-1).*vn(2:N-1);
                for i = 2:N-1
                    DF(s_vn+i-1,s_v0) = 2*dx*eta*thetan(i);
                    DF(s_vn+i-1,s_vn+i-1) = -2*dx*thetan(i)*(eta+nust(i));
                    DF(s_vn+i-1,s_tn+i-2) = ksigman;
                    DF(s_vn+i-1,s_tn+i-1) = 2*dx*eta*(v0-vn(i)) - 2*dx*nust(i)*vn(i);
                    DF(s_vn+i-1,s_tn+i)   = -ksigman;
                end
                
                Fn(s_vn+N-1) = thetan(N)*(vn(N) - v0) + Jactinf;
                DF(s_vn+N-1,s_vn+N-1) = thetan(N);
                DF(s_vn+N-1,s_tn+N-1) = vn(N) - v0;
                DF(s_vn+N-1,s_tc+N-1) = DJactinfDthcN;
                DF(s_vn+N-1,s_v0) = -thetan(N);
                
                %% Equations for thetan and the Derivatives
                Fn(s_tn) = thetan(1)*(vn(1)-v0);
                DF(s_tn,s_vn) = thetan(1);
                DF(s_tn,s_tn) = vn(1)-v0;
                DF(s_tn,s_v0) = -thetan(1);
                
                Fn(s_tn+1:s_tn+N-2) = (thetan(3:N).*vn(3:N)-thetan(1:N-2).*vn(1:N-2))...
                    - v0*(thetan(3:N)-thetan(1:N-2)) + 2*dx*gamma(2:N-1).*thetan(2:N-1);
                for i = 2:N-1
                    DF(s_tn+i-1,[s_vn+i-2,s_vn+i]) = [-thetan(i-1), thetan(i+1)];
                    DF(s_tn+i-1,s_tn+i-2) = -vn(i-1) + v0;
                    DF(s_tn+i-1,s_tn+i-1) = 2*dx*gamma(i);
                    DF(s_tn+i-1,s_tn+i)   =  vn(i+1) - v0;
                    DF(s_tn+i-1,s_v0) = -(thetan(i+1) - thetan(i-1));
                end
                
                Fn(s_tn+N-1) = (thetan(N)*vn(N)-thetan(N-1)*vn(N-1)) ...
                    - v0*(thetan(N)-thetan(N-1)) + dx*gamma(N)*thetan(N);
                DF(s_tn+N-1,[s_vn+N-2,s_vn+N-1]) = [-thetan(N-1), thetan(N)];
                DF(s_tn+N-1,s_tn+N-2) = -vn(N-1) + v0;
                DF(s_tn+N-1,s_tn+N-1) =  vn(N)   - v0 + dx*gamma(N);
                DF(s_tn+N-1,s_v0) = -(thetan(N)-thetan(N-1));
                
                %% Equations for thetac and the Derivatives
                Fn(s_tc) = thetac(2) - thetac(1);
                DF(s_tc,s_tc)   = -1;
                DF(s_tc,s_tc+1) =  1;
                
                Fn(s_tc+1:s_tc+N-2) = thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N)...
                    + dx^2/Dtc*gamma(2:N-1).*thetan(2:N-1);
                for i = 2:N-1
                    DF(s_tc+i-1,s_tn+i-1) = dx^2/Dtc*gamma(i);
                    DF(s_tc+i-1,s_tc+i-2) = 1;
                    DF(s_tc+i-1,s_tc+i-1) = -2;
                    DF(s_tc+i-1,s_tc+i)   = 1;
                end
                
                Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) + sum(thetan(2:N)+thetac(2:N)))...
                    - L*(Thetan + Thetac)/dx;
                DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
                DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
                DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
                DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;
                
                %% Equation for v0 and the Derivatives
                if FadMode == 1
                    Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad)*v0 ...
                        + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                    DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                    DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                    DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                    DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
                    DF(s_v0,s_v0) = dg + kad;
                elseif FadMode == 2
                    Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad*thetan(1))*v0 ...
                        + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                    DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                    DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                    DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                    DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
                    DF(s_v0,s_tn+1) = kad*v0;
                    DF(s_v0,s_v0) = dg + kad*thetan(1);
                end
                
                temp_Fn = Fn;
                
                %% Solve for the matrix
                DF = sparse(DF);
                X = X - DF\Fn;
                
                if sum(isnan(X)) >= 1 || iter == Iter || norm(temp_Fn) < norm(Fn)
                    Solution = 0;
                    vn     = NaN(N,1);
                    thetan = NaN(N,1);
                    thetac = NaN(N,1);
                    mn     = NaN(N,1);
                    mc     = NaN(N,1);
                    v0     = NaN;
                    break
                else
                    Solution = 1;
                end
                
                vn = X(s_vn:s_vn+N-1);
                thetan = X(s_tn:s_tn+N-1);
                thetac = X(s_tc:s_tc+N-1);
                v0 = X(s_v0);
                
                if iter > 1
                    error = abs((X-temp_X)./X);
                    error = sum(error)/(N_var);
                    if error < 1d-5 || iter == Iter
                        ITER = false;
                    end
                end
                temp_X = X;
            end
            
            if loop1 == 1 && ig == 1
                temp_vn_loop1 = vn;
                temp_thetan_loop1 = thetan;
                temp_thetac_loop1 = thetac;
                temp_v0_loop1 = v0;
            end
            if ig == 1
                temp_vn_ig = vn;
                temp_thetan_ig = thetan;
                temp_thetac_ig = thetac;
                temp_v0_ig = v0;
            end
            
%             lambda = sqrt((eta+etast)/ksigman)*gamma0.^(1/2);
%             int_vn = -1/2*sum(vn(1:N-1)+vn(2:N))*dx;
%             tauD = (eta+etast)/ksigman/lambda*L;
%             
%             dthetan = [thetan(2)-thetan(1);(thetan(3:N)-thetan(1:N-2))/2;thetan(N)-thetan(N-1)]/dx;
            
            V0(ig) = v0;
            
        end
        [~,ig] = max(V0);
        MAXGAMMA(loop1,loop2) = Gamma(ig);
        
        
        
        if loop1 == N1 && loop2 == N2
            figure(21);
            subplot(3,1,1)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,gamma,x*1d-3,ksigman*thetan);
            xlabel('x ({\mu}m)','fontsize',15)
            ylabel(hAx(1),'\gamma (1/s)','fontsize',15) % left y-axis
            ylabel(hAx(2),'\sigma (Pa)','fontsize',15) % right y-axis
            set(gca,'fontsize',15); set(hAx(2),'fontsize',15);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(3,1,2)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,thetan*1d3,x*1d-3,thetac*1d3);
            xlabel('x ({\mu}m)','fontsize',15)
            ylabel(hAx(1),'\theta_n ({\mu}M)','fontsize',15) % left y-axis
            ylabel(hAx(2),'\theta_c ({\mu}M)','fontsize',15) % right y-axis
            set(gca,'fontsize',15); set(hAx(2),'fontsize',15);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(3,1,3)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,vn,x*1d-3,vn-v0);
            xlabel('x ({\mu}m)','fontsize',15)
            ylabel(hAx(1),'v_n (nm/s)','fontsize',15) % left y-axis
            ylabel(hAx(2),'v_n - v_0 (nm/s)','fontsize',15) % right y-axis
            set(gca,'fontsize',15); set(hAx(2),'fontsize',15);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
        end
        
        
    end
end

%% Plotting

if Parameter1 == 0 && Parameter2 == 0
    
    figure(103)
    plot(Gamma,V0,'-','linewidth',2); hold on
    set(gca,'fontsize',15);
    xlabel('\gamma (1/s)','fontsize',15)
    ylabel('v_0 (nm/s)','fontsize',15)
    xlim([0 L*1d-3])
    box off
%     legend('\gamma = 5\times10^{-4} 1/s','\gamma = 1\times10^{-3} 1/s',...
%         '\gamma = 5\times10^{-3} 1/s','\gamma = 1\times10^{-2} 1/s')
    
    
end

if Parameter1 == 0 && Parameter2 > 0
    figure(201)
    plot(Xm,MAXGAMMA,'-','linewidth',2); hold on
    AxesProperties20201003
    ylabel('\gamma (1/s)')
    box off
    %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
    %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
    grid on
    
end


if Parameter1 > 0 && Parameter2 > 0
    [Xmesh,Ymesh] = meshgrid(Xm,Ym);
    data = log10(MAXGAMMA);
    figure(301)
    [~,hd] = contour(Xmesh,Ymesh,data,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(data)),max(max(data)),25));
    title('log_{10}\gamma (1/s)','fontsize',15,'FontWeight','normal')
    AxesProperties20201003
    
end




