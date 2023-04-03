% This is part of the orginal codes used in the following paper:
% http://www.molbiolcell.org/cgi/doi/10.1091/mbc.E22-10-0494
% On the role of myosin-induced actin depolymerization during cell migration
% If you have questions, feel free to contact Dr. Yizeng Li.


% Two-phase: cytosol and actomyosin network
% Steady-state
% Unknowns:
% vn, theta_n, theta_c, m_n, m_c, v0
% The velocities are written in a fixed frame but the equations are written
% in the frame of the cell.


clear
clc
%close all

Parameter1 = 3;
Parameter2 = 9;
% 1: Jactinf
% 2: nust
% 3: gamma0
% 4: kad
% 5: dg
% 6: eta
% 7: Dtc
% 8: ksigman
% 9: ksigmaa
% 10: fextf
% 11: koff
% 12: kon

FadMode = 1;
% FadMode = 1: Fad = kad*v0
% FadMode = 2: Fad = kad*thetan(1)*v0

JactinMode = 2;
% JactinMode = 1: Jactinf = Jactinf0*thetac(N)
% JactinMode = 2: Jactinf = Jactinf0*thetac(N)/(thetacc + thetac(N))

GammaMode = 1;
% GammaMode = 1: gamma = constant
% GammaMode = 2: gamma = gamma0 + gammaa*sigma_a
% GammaMode = 3: gamma = gamma0 + gammaa*sigma_a/(sigma_ac + sigma_a)



%% Parameters are in units: nm, s, Pa & mM

L = 50.d3;              % (nm) cell length
p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell

Thetac = 0.1;           % (mM) reference value of G-actin
Thetan = 0.2;           % (mM) reference value of F-actin
Theta  = Thetac + Thetan; % (mM) reference total G- and F-actin, \int_0^L (thetan + thetac)dx/L

Mc = 5d-3;              % (mM) reference value of free myosin
Mn = 5d-3;              % (mM) reference value of

if JactinMode == 1
    Jactinf0 = 50;          % (nm/s) Jactin^f = Jactinf0*thetac^f
elseif JactinMode == 2
    Jactinf0 = 6;           % (nm mM/s) Jactinf = Jactinf0*thetac^f/(thetacc + thetac^f)
    thetacc  = 0.2d-3;         % (mM) Critical value for actin polymerization
end

if GammaMode == 1
    gamma0 = 5d-4;          % (1/s) constant rate of actin depolymerization
    %gamma0 = 3d-3;          % (1/s) constant rate of actin depolymerization
    %gamma0 = 3d-3;
    % gamma0 = 1d-3;
    % gamma0 = 5d-3;
    % gamma0 = 1d-2;
elseif GammaMode == 2
    gamma0 = 1d-4;      % (1/s) gamma = gamma0 + gammmaa*sigma_a
    gammaa = 2d-3;      % (1/Pa/s)
elseif GammaMode == 3
    gamma0 = 5d-4;      % (1/s) gamma = gamma0 + gammaa*sigma_a/(sigma_ac + sigma_a)
    gammaa = 4d-3;      % (1/s)
    sigma_ac = 1;       % (Pa)
end

kon = 1d-0;     % (1/s mM) ka*m_c*theta_n
koff = 1.5d-0;   % (1/s) ko*m_n 

ksigman = 100*1d2;          % (Pa /mM) Coefficient of passive actin network pressure
ksigmaa = 5d2;          % (Pa /mM) Coefficient of active actin network contraction
%ksigmaa = 1d6; 

etast0 = 100*1d-4;           % (Pa s/nm^2/mM)
eta   = 1d-8;           % (Pa s/nm^2/mM)
dg    = 0d-6;           % (Pa s/nm) coefficient of hydraulic resistance
if FadMode == 1
    kad = 100*3d-1;           % (Pa s/nm) adhesive force, Fad^b = kad*v0
elseif FadMode == 2
    kad = 1d0;          % (Pa s/mM/nm) adhesive force, Fad^b = kad*thetan^b*v0
end

Dtc = 1.d7;             % (nm^2/s) diffusion constant for theta_c
Dmc = 1.d6;             % (nm^2/s) diffusion constant for m_c
Dmn = 1.d5;             % (nm^2/s) diffusion constant for m_n

%%
Iter = 21;
N  = 101;
dx = L/(N-1);
x  = linspace(0,L,N);

% number of variables of the model
n_var = 6; % v_n, theta_n, theta_c, m_n, m_c, v0
var_legh = [N,N,N,N,N,1]; % length of eavh variable
var_strt = zeros(1,n_var); % start position of each variable
var_strt(1) = 1;
for in = 2:n_var
    var_strt(in) = var_strt(in-1) + var_legh(in-1);
end
N_var = var_strt(n_var) + var_legh(n_var) - 1; % total size of the matrix

s_vn = var_strt(1); l_vn = var_legh(1);
s_tn = var_strt(2); l_tn = var_legh(2);
s_tc = var_strt(3); l_tc = var_legh(3);
s_mn = var_strt(4); l_mn = var_legh(4);
s_mc = var_strt(5); l_mc = var_legh(5);
s_v0 = var_strt(6); l_v0 = var_legh(6);


Selection = 1;
Sweep_full_20230113;

V0  = zeros(N1,N2);
VNF = zeros(N1,N2);
TNB = zeros(N1,N2);
TNF = zeros(N1,N2);
FAD = zeros(N1,N2);
POWER = zeros(N1,N2);
DIFFERENCE = zeros(N1,N2);

Solution = 0;
for loop2 = 1:N2
    loop2
    Selection = 2;
    Sweep_full_20230113;
    for loop1 = 1:N1
        Selection = 3;
        Sweep_full_20230113;
        
        etast = etast0*logspace(-0,0,N)';
        
        sigma_a = ksigmaa*Mn;
        if GammaMode == 1
            gamma = gamma0;
        elseif GammaMode == 2
            gamma = gamma0 + gammaa*sigma_a;
        elseif GammaMode == 3
            gamma = gamma0 + gammaa*sigma_a./(sigma_ac + sigma_a);
        end
        

        % initial guess
        if (loop1 == 1 && loop2 == 1) || Solution == 0
%             thetan = linspace(Thetan*1,Thetan*1,N)';
%             thetac = linspace(Thetac*1,Thetac*1,N)';
%             mn = Mn*ones(N,1);
%             mc = Mc*ones(N,1);
%             v0 = ((1/2)*(etast0*Theta*L)*(Jactinf0*gamma*L)...
%                 -fextf*(Jactinf0+gamma*L))./((etast0*Theta*L)*Jactinf0 ...
%                 + (kad*Thetac+dg)*(Jactinf0+gamma*L));
%             vn = linspace(v0,v0-gamma*L,N)';

            load('Initial_myosin_full_ksigmaa1d3_20230116.mat');
        elseif loop1 == 1 && loop2 > 1
            vn = temp_vn;
            thetan = temp_thetan;
            thetac = temp_thetac;
            mn = temp_mn;
            mc = temp_mc;
            v0 = temp_v0;
        end
        X = [vn; thetan; thetac; mn; mc; v0];
        
        iter = 0;
        ITER = true;
        while ITER
            iter = iter + 1;
            
            DF = zeros(N_var,N_var);
            Fn = zeros(N_var,1);
            
            sigma_n = ksigman*thetan;     % N by 1 vector
            sigma_a = ksigmaa*mn;         % N by 1 vecyor
            sigma   = sigma_n - sigma_a;
            dsigmadtn =  ksigman;
            dsigmadmn = -ksigmaa;
            
            if GammaMode == 1
                gamma = gamma0*logspace(-0,0,N)';
                dgammadtn = zeros(N,1);
                dgammadmn = zeros(N,1);
            elseif GammaMode == 2
                gamma = gamma0 + gammaa*sigma_a;
                dgammadtn = zeros(N,1);
                dgammadmn = gammaa*ksigmaa*logspace(-0,0,N)';
            elseif GammaMode == 3
                gamma = gamma0 + gammaa*sigma_a./(sigma_ac + sigma_a);
                dgammadtn = zeros(N,1);
                dgammadmn = gammaa*ksigmaa*sigma_ac./(sigma_ac + sigma_a).^2;
            end
            
            if JactinMode == 1
                Jactinf = Jactinf0*thetac(N);
                DJactinfDtcN = Jactinf0;
            elseif JactinMode == 2
                Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
                DJactinfDtcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
            end
            
            %% Equations for vn and the Derivatives
            
            Fn(s_vn) = thetan(N)*(vn(N) - v0) + Jactinf;
            DF(s_vn,s_vn+N-1) = thetan(N);
            DF(s_vn,s_tn+N-1) = vn(N) - v0;
            DF(s_vn,s_tc+N-1) = DJactinfDtcN;
            DF(s_vn,s_v0) = -thetan(N);
            
            
            Fn(s_vn+1:s_vn+N-1) = -(sigma(2:N)-sigma(1:N-1)) + dx*eta*thetan(2:N).*(v0-vn(2:N)) ...
                - dx*etast(2:N).*thetan(2:N).*vn(2:N);
            for i = 2:N
                DF(s_vn+i-1,s_v0) = dx*eta*thetan(i);
                DF(s_vn+i-1,s_vn+i-1) = -dx*thetan(i)*(eta+etast(i));
                DF(s_vn+i-1,s_tn+i-2) = dsigmadtn;
                DF(s_vn+i-1,s_tn+i-1) = -dsigmadtn + dx*eta*(v0-vn(i)) - dx*etast(i)*vn(i);
                DF(s_vn+i-1,[s_mn+i-2,s_mn+i-1]) = [1,-1]*dsigmadmn;
            end
            
            
            %% Equations for thetan and the Derivatives
            Fn(s_tn) = thetan(1)*(vn(1)-v0);
            DF(s_tn,s_vn) = thetan(1);
            DF(s_tn,s_tn) = vn(1)-v0;
            DF(s_tn,s_v0) = -thetan(1);
            
            Fn(s_tn+1:s_tn+N-1) = (thetan(2:N).*vn(2:N)-thetan(1:N-1).*vn(1:N-1))...
                - v0*(thetan(2:N)-thetan(1:N-1)) + dx*gamma(2:N).*thetan(2:N);
            for i = 2:N
                DF(s_tn+i-1,[s_vn+i-2,s_vn+i-1]) = [-thetan(i-1), thetan(i)];
                DF(s_tn+i-1,s_tn+i-2) = -vn(i-1) + v0;
                DF(s_tn+i-1,s_tn+i-1) = dx*gamma(i) + dx*thetan(i)*dgammadtn(i) ...
                    + vn(i) - v0;
                DF(s_tn+i-1,s_mn+i-1) = dx*thetan(i)*dgammadmn(i);
                DF(s_tn+i-1,s_v0) = -(thetan(i) - thetan(i-1));
            end
            
            
            %% Equations for thetac and the Derivatives
            Fn(s_tc) = thetac(2) - thetac(1);
            DF(s_tc,s_tc)   = -1;
            DF(s_tc,s_tc+1) =  1;
            
            Fn(s_tc+1:s_tc+N-2) = thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N)...
                + dx^2/Dtc*gamma(2:N-1).*thetan(2:N-1);
            for i = 2:N-1
                DF(s_tc+i-1,s_tn+i-1) = dx^2/Dtc*gamma(i) + dx^2/Dtc*thetan(i)*dgammadtn(i);
                DF(s_tc+i-1,s_tc+i-2) = 1;
                DF(s_tc+i-1,s_tc+i-1) = -2;
                DF(s_tc+i-1,s_tc+i)   = 1;
                DF(s_tc+i-1,s_mn+i-1) = dx^2/Dtc*thetan(i)*dgammadmn(i);
            end
            
            Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) + sum(thetan(2:N)+thetac(2:N)))...
                - L*(Thetan + Thetac)/dx;
            DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
            DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
            DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
            DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;
            
            %% Equations for mn and the Derivatives
            Fn(s_mn) = mn(1)*(vn(1)-v0) - Dmn/dx*(mn(2)-mn(1));
            DF(s_mn,s_vn)   = mn(1);
            DF(s_mn,s_mn)   = (vn(1) - v0) + Dmn/dx;
            DF(s_mn,s_mn+1) = - Dmn/dx;
            DF(s_mn,s_v0)   = -mn(1);
            
            Fn(s_mn+1:s_mn+N-2) = (vn(3:N).*mn(3:N) - vn(1:N-2).*mn(1:N-2)) ...
                - v0*(mn(3:N) - mn(1:N-2)) ...
                - 2*Dmn/dx*(mn(1:N-2) - 2*mn(2:N-1) + mn(3:N)) ...
                - 2*dx*kon*mc(2:N-1).*thetan(2:N-1) + 2*dx*koff*mn(2:N-1);
            for i = 2:N-1
                DF(s_mn+i-1,[s_vn+i-2,s_vn+i]) = [-mn(i-1), mn(i+1)];
                DF(s_mn+i-1,s_tn+i-1) = -2*dx*kon*mc(i);
                DF(s_mn+i-1,s_mn+i-2) = -vn(i-1) + v0 - 2*Dmn/dx;
                DF(s_mn+i-1,s_mn+i-1) = 4*Dmn/dx + 2*dx*koff;
                DF(s_mn+i-1,s_mn+i)   =  vn(i+1) - v0 - 2*Dmn/dx;
                DF(s_mn+i-1,s_mc+i-1) = -2*dx*kon*thetan(i);
                DF(s_mn+i-1,s_v0)     = -(mn(i+1) - mn(i-1));
            end
            
            Fn(s_mn+N-1) = mn(N)*(vn(N)-v0) - Dmn/dx*(mn(N)-mn(N-1));
            DF(s_mn+N-1,s_vn+N-1) = mn(N);
            DF(s_mn+N-1,s_mn+N-2) = Dmn/dx;
            DF(s_mn+N-1,s_mn+N-1) = (vn(N) - v0) - Dmn/dx;
            DF(s_mn+N-1,s_v0)     = -mn(N);
            
            %% Equations for mc and the Derivatives
            Fn(s_mc) = mc(2)-mc(1);
            DF(s_mc,s_mc)   = -1;
            DF(s_mc,s_mc+1) = 1;
            
            Fn(s_mc+1:s_mc+N-2) = -2*Dmc/dx*(mc(1:N-2) - 2*mc(2:N-1) + mc(3:N)) ...
                + 2*dx*kon*mc(2:N-1).*thetan(2:N-1) - 2*dx*koff*mn(2:N-1);
            for i = 2:N-1
                DF(s_mc+i-1,s_tn+i-1) = 2*dx*kon*mc(i);
                DF(s_mc+i-1,s_mn+i-1) = -2*dx*koff;
                DF(s_mc+i-1,s_mc+i-2) = - 2*Dmc/dx;
                DF(s_mc+i-1,s_mc+i-1) = 4*Dmc/dx + 2*dx*kon*thetan(i);
                DF(s_mc+i-1,s_mc+i)   = - 2*Dmc/dx;
            end
            
            Fn(s_mc+N-1) = 1/2*(sum(mn(1:N-1)+mc(1:N-1)) + sum(mn(2:N)+mc(2:N))) - L*(Mn + Mc)/dx;
            DF(s_mc+N-1,[s_mn,s_mn+N-1]) = 1/2;
            DF(s_mc+N-1,s_mn+1:s_mn+N-2) = 1;
            DF(s_mc+N-1,[s_mc,s_mc+N-1]) = 1/2;
            DF(s_mc+N-1,s_mc+1:s_mc+N-2) = 1;
            
            %% Equation for v0 and the Derivatives
            if FadMode == 1
                Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad)*v0 ...
                    + dx/2*(sum(etast(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
                    + sum(etast(2:N).*thetan(2:N).*vn(2:N)));
                DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[etast(1)*thetan(1),etast(N)*thetan(N)];
                DF(s_v0,s_vn+1:s_vn+N-2) = dx*etast(2:N-1).*thetan(2:N-1);
                DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[etast(1)*vn(1),etast(N)*vn(N)];
                DF(s_v0,s_tn+1:s_tn+N-2) = dx*etast(2:N-1).*vn(2:N-1);
                DF(s_v0,s_v0) = dg + kad;
            elseif FadMode == 2
                Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad*thetan(1))*v0 ...
                    + dx/2*(sum(etast(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
                    + sum(etast(2:N).*thetan(2:N).*vn(2:N)));
                DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[etast(1)*thetan(1),etast(N)*thetan(N)];
                DF(s_v0,s_vn+1:s_vn+N-2) = dx*etast(2:N-1).*thetan(2:N-1);
                DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[etast(1)*vn(1),etast(N)*vn(N)];
                DF(s_v0,s_tn+1:s_tn+N-2) = dx*etast(2:N-1).*vn(2:N-1);
                DF(s_v0,s_tn) = DF(s_v0,s_tn) + kad*v0;
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
            mn = X(s_mn:s_mn+N-1);
            mc = X(s_mc:s_mc+N-1);
            v0 = X(s_v0);
            
            if iter > 1
                error = abs((X-temp_X)./X);
                error = sum(error)/(N_var);
                if error < 1d-6 || iter == Iter
                    ITER = false;
                end
            end
            temp_X = X;
        end
        
        if loop1 == 1
            temp_vn = vn;
            temp_thetan = thetan;
            temp_thetac = thetac;
            temp_mn = mn;
            temp_mc = mc;
            temp_v0 = v0;
        end
        
        V0(loop1,loop2) = v0;
        VNF(loop1,loop2) = vn(N);
        TNB(loop1,loop2) = thetan(1)*1d3;
        TNF(loop1,loop2) = thetan(N)*1d3;
        FAD(loop1,loop2) = kad*thetac(1)*v0;
        POWER(loop1,loop2) = v0*fextf;
        DIFFERENCE(loop1,loop2) = (Jactinf ...
            - dx/2*(sum(gamma(1:N-1).*thetan(1:N-1))+sum(gamma(2:N).*thetan(2:N))))/Jactinf*100;
        
        if loop1 == N1 && loop2 == N2 
            figure(22);
            subplot(4,1,1)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,gamma,x*1d-3,ksigman*thetan);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'\gamma (1/s)','fontsize',18) % left y-axis
            ylabel(hAx(2),'\sigma_n (Pa)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(4,1,2)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,thetan*1d3,x*1d-3,thetac*1d3);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'\theta_n ({\mu}M)','fontsize',18) % left y-axis
            ylabel(hAx(2),'\theta_c ({\mu}M)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(4,1,3)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,mn*1d3,x*1d-3,mc*1d3);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'m_n ({\mu}M)','fontsize',18) % left y-axis
            ylabel(hAx(2),'m_c ({\mu}M)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(4,1,4)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,vn,x*1d-3,vn-v0);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'v_n (nm/s)','fontsize',18) % left y-axis
            ylabel(hAx(2),'v_n - v_0 (nm/s)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);

        end
        
    end
end

%% Plotting

if Parameter1 == 0 && Parameter2 == 0
    fprintf('v0 = %4.4f nm/s\n',v0);
    fprintf('Jactinf = %4.4f nm/s\n',Jactinf0*thetac(N)/thetan(N));
    
    figure(102)
    plot(x*1d-3,thetan*1d3,'-','linewidth',2); hold on
    set(gca,'fontsize',18); 
    xlabel('x - x^b({\mu}m)','fontsize',18)
    ylabel('\theta_n ({\mu}M)','fontsize',18)
    xlim([0 L*1d-3])
    box off
    % legend('\gamma = 5\times10^{-4} 1/s','\gamma = 1\times10^{-3} 1/s',...
    %     '\gamma = 5\times10^{-3} 1/s','\gamma = 1\times10^{-2} 1/s')
    
    % legend('k_{\sigma_a} = 1\times10^{2} Pa/mM','k_{\sigma_a} = 5\times10^{2} Pa/mM',...
    %     'k_{\sigma_a} = 1\times10^{3} Pa/mM','k_{\sigma_a} = 1\times10^{3} Pa/mM');
    
%     legend('k_{\sigma_a} = 5\times10^{1} Pa/mM','k_{\sigma_a} = 1\times10^{2} Pa/mM',...
%        'k_{\sigma_a} = 5\times10^{2} Pa/mM','k_{\sigma_a} = 1\times10^{3} Pa/mM');
    
    
    figure(103)
    plot(x*1d-3,mn*1d3,'-','linewidth',2); hold on
    set(gca,'fontsize',18); 
    xlabel('x - x^b ({\mu}m)','fontsize',18)
    ylabel('m_n ({\mu}M)','fontsize',18)
    xlim([0 L*1d-3])
    box off
    % legend('\gamma = 5\times10^{-4} 1/s','\gamma = 1\times10^{-3} 1/s',...
    %    '\gamma = 5\times10^{-3}1/s','\gamma = 1\times10^{-2} 1/s')
%     legend('k_{\sigma_a} = 2\times10^{2} Pa/mM','k_{\sigma_a} = 1\times10^{3} Pa/mM',...
%         'k_{\sigma_a} = 5\times10^{3} Pa/mM','k_{\sigma_a} = 1\times10^{4} Pa/mM');
%     legend('k_{\sigma_a} = 1\times10^{2} Pa/mM','k_{\sigma_a} = 1\times10^{3} Pa/mM',...
%         'k_{\sigma_a} = 1\times10^{4} Pa/mM','k_{\sigma_a} = 1\times10^{5} Pa/mM');

    figure(104)
    plot(x*1d-3,vn,'-','linewidth',2); hold on
    set(gca,'fontsize',18);
    xlabel('x - x^b ({\mu}m)','fontsize',18)
    ylabel('v_n (nm/s)','fontsize',18) 
    xlim([0 L*1d-3])
    box off
%     legend('k_{\sigma_a} = 2\times10^{2} Pa/mM','k_{\sigma_a} = 1\times10^{3} Pa/mM',...
%         'k_{\sigma_a} = 5\times10^{3} Pa/mM','k_{\sigma_a} = 1\times10^{4} Pa/mM');
    
end

if Parameter1 == 0 && Parameter2 > 0
    figure(201)
    index = V0 >= 0;
    plot(Xm(index),V0(index),'-','linewidth',2); hold on
    AxesProperties20201009
    ylabel('v_0 (nm/s)')
    box off
%     legend('k_{\sigma_n} = 1\times10^4 Pa/mM','k_{\sigma_n} = 2\times10^4 Pa/mM',...
%         'k_{\sigma_n} = 3\times10^4 Pa/mM','k_{\sigma_n} = 4\times10^4 Pa/mM')
    
    % legend('k_{on} = 0.5 /s mM','k_{on} = 1 /s mM','k_{on} = 2 /s mM','k_{on} = 5 /s mM');
    
    
    if Parameter2 == 10 
        figure(202)
        index = V0 >= 0;
        plot(V0(index),POWER(index)*1d-3,'-','linewidth',2); hold on
        set(gca,'fontsize',18);
        xlabel('v_0 (nm/s)')
        ylabel('Power ({\mu}W/m^2)')
        box off
        % legend('k_{\sigma_a} = 2\times10^{3} Pa/mM','k_{\sigma_a} = 7\times10^{3} Pa/mM');
    end
    
end

if Parameter1 > 0 && Parameter2 > 0
    [Xmesh,Ymesh] = meshgrid(Xm,Ym);
    
    figure(301)
    [~,hd] = contour(Xmesh,Ymesh,V0,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(V0)),max(max(V0)),31));
    title('v_0 (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    
    figure(302)
    subplot(2,2,1)
    [~,hd] = contour(Xmesh,Ymesh,V0,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(V0)),max(max(V0)),25));
    title('v_0 (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    subplot(2,2,2)
    [~,hd] = contour(Xmesh,Ymesh,VNF,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(VNF)),max(max(VNF)),25));
    title('v_n^f (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    subplot(2,2,3)
    [~,hd] = contour(Xmesh,Ymesh,TNB,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(TNB)),max(max(TNB)),25));
    title('\theta_n^b ({\mu}M)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    subplot(2,2,4)
    [~,hd] = contour(Xmesh,Ymesh,TNF,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(TNF)),max(max(TNF)),25));
    title('\theta_n^f ({\mu}M)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    
    if Parameter2 == 10 
        figure(303)
        index = POWER < -inf;
        POWER(index) = NaN;
        [~,hd] = contour(Xmesh,Ymesh,POWER,'fill','on'); colorbar;
        set(hd,'LevelList',linspace(min(min(POWER)),max(max(POWER)),25));
        title('Power (nW/m^2)','fontsize',18,'FontWeight','normal')
        AxesProperties20201009
        
%         TargetKsigmaa = [1d2, 3d2,5d2,2d3];
%         nT = length(TargetKsigmaa);
%         for inT = 1:nT
%             [~,idx] = min(abs(Ym-TargetKsigmaa(inT)));
%             
%             figure(3031)
%             vcell = V0(idx,:);
%             index = vcell >=0;
%             plot(Xm(index),vcell(index),'-','linewidth',2); hold on
%             set(gca,'fontsize',18);
%             xlabel('f_{stall} (Pa/m^2)')
%             ylabel('v_0 (nm/s)')
%             box off
%             legend('k_{\sigma_a} = 1\times10^{2} Pa/mM','k_{\sigma_a} = 3\times10^{2} Pa/mM',...
%                 'k_{\sigma_a} = 5\times10^{2} Pa/mM','k_{\sigma_a} = 2\times10^{3} Pa/mM');
%         end
        
        TargetKsigmaa = [1d3, 1.5d4];
        nT = length(TargetKsigmaa);
        for inT = 1:nT
            [~,idx] = min(abs(Ym-TargetKsigmaa(inT)));
            figure(3033)
            vcell = V0(idx,:);
            power = POWER(idx,:);
            index = vcell >=0;
            plot(vcell(index),power(index),'-','linewidth',2); hold on
            set(gca,'fontsize',18);
            xlabel('v_0 (nm/s)')
            ylabel('Power (nW/m^2)')
            box off
            legend('k_{\sigma_a} = 1\times10^{3} Pa/mM','k_{\sigma_a} = 1.5\times10^{4} Pa/mM');
        end
    end
    
    if Parameter1 == 10
        
        figure(304)
        index = POWER < -inf;
        POWER(index) = NaN;
        [~,hd] = contour(Xmesh,Ymesh,POWER,'fill','on'); colorbar;
        set(hd,'LevelList',linspace(min(min(POWER)),max(max(POWER)),25));
        title('Power (nW/m^2)','fontsize',18,'FontWeight','normal')
        AxesProperties20201009
        
        Power = max(POWER);
        figure(305)
        plot(Xm,Power*1d-3,'-','linewidth',2); hold on
        AxesProperties20201009
        ylabel('Max. Power ({\mu}W/m^2)')
        vaxis = axis;
        vaxis(3) = min(Power)*1d-3;
        vaxis(4) = max(Power)*1d-3;
        axis(vaxis);
        box off
    end
    
    figure(306)
    [~,hd] = contour(Xmesh,Ymesh,DIFFERENCE,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(DIFFERENCE)),max(max(DIFFERENCE)),25));
    title('(J_{actin}^f - \int_0^L{\gamma}\theta_ndx)/J_{actin}^f (%)',...
        'fontsize',18,'FontWeight','normal')
    AxesProperties20201009
    
    figure(307)
    [~,hd] = contour(Xmesh,Ymesh,FAD,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(FAD)),max(max(FAD)),25));
    title('F_{ad} (Pa)','fontsize',18,'FontWeight','normal')
    AxesProperties20201009

end
