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
% 3: gamma
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

Parameter1 = 0;
Parameter2 = 0;

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
    Jactinf0 = 30;          % (nm/s) coefficient of actin polymerization
elseif JactinMode == 2
    Jactinf0 = 6;           % (nm mM/s) coefficient of actin polymerization
    thetacc  = 0.2d-3;         % (mM) Critical value for actin polymerization
end
 gamma0 = 5d-4;          % (1/s) 
  gamma0 = 1d-3;
  gamma0 = 5d-3;
  gamma0 = 1d-2;

ksigman = 100*1d2;          % (Pa /mM) Coefficient of passive actin network pressure

etast = 100*1d-4;           % (Pa s/nm^2/mM)
eta   = 1d-8;           % (Pa s/nm^2/mM)
dg    = 0*1d-6;           % (Pa s/nm) coefficient of hydraulic resistance
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

VNF = zeros(N1,N2);
VC = zeros(N1,N2);
V0 = zeros(N1,N2);
TNB = zeros(N1,N2);
TNF = zeros(N1,N2);
JWATERF = zeros(N1,N2);
TAUF = zeros(N1,N2);
TAUB = zeros(N1,N2);
PSTARF = zeros(N1,N2);
PSTARB = zeros(N1,N2);
PCF = zeros(N1,N2);
PCB = zeros(N1,N2);
V0_L = zeros(N1,N2);
POWER = zeros(N1,N2);
VARIATION1 = zeros(N1,N2);
VARIATION2 = zeros(N1,N2);
VARIATION3 = zeros(N1,N2);
VARIATION4 = zeros(N1,N2);
VARIATION5 = zeros(N1,N2);
DIFFERENCE = zeros(N1,N2);

Solution = 0;
for loop2 = 1:N2
    loop2
    Selection = 2;
    Sweep_20201003;
    for loop1 = 1:N1
        Selection = 3;
        Sweep_20201003;
        
        gamma = gamma0*logspace(-0,0,N)';
        nust = etast*logspace(-0,0,N)';

        % initial guess
        if (loop1 == 1 && loop2 == 1) || Solution == 0
            thetan = linspace(Theta*1,Theta*1,N)'*(Jactinf0/(Jactinf0+gamma0*L));
            thetac = linspace(Theta*1,Theta*1,N)'*(gamma0*L/(Jactinf0+gamma0*L));
            v0 = (1/2)*(etast*Theta*L)*(Jactinf0*gamma0*L)./((etast*Theta*L)*Jactinf0 ...
                + (kad+dg)*(Jactinf0+gamma0*L));
            vn = linspace(v0,v0-gamma0*L,N)';
        elseif loop1 == 1 && loop2 > 1
            vn = temp_vn;
            thetan = temp_thetan;
            thetac = temp_thetac;
            v0 = temp_v0;
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
                    + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
                    + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
                DF(s_v0,s_v0) = dg + kad;
            elseif FadMode == 2
                Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dg + kad*thetan(1))*v0 ...
                    + dx/2*(sum(nust(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
                    + sum(nust(2:N).*thetan(2:N).*vn(2:N)));
                DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[nust(1)*thetan(1),nust(N)*thetan(N)];
                DF(s_v0,s_vn+1:s_vn+N-2) = dx*nust(2:N-1).*thetan(2:N-1);
                DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[nust(1)*vn(1),nust(N)*vn(N)];
                DF(s_v0,s_tn+1:s_tn+N-2) = dx*nust(2:N-1).*vn(2:N-1);
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
        
        if loop1 == 1
            temp_vn = vn;
            temp_thetan = thetan;
            temp_thetac = thetac;
            temp_v0 = v0;
        end
        
        lambda = sqrt((eta+etast)/ksigman)*gamma0.^(1/2);
        int_vn = -1/2*sum(vn(1:N-1)+vn(2:N))*dx;
        tauD = (eta+etast)/ksigman/lambda*L; 
        
        dthetan = [thetan(2)-thetan(1);(thetan(3:N)-thetan(1:N-2))/2;thetan(N)-thetan(N-1)]/dx;
 
        
        V0(loop1,loop2) = v0;
        VNF(loop1,loop2) = vn(N);
        TNB(loop1,loop2) = thetan(1);
        TNF(loop1,loop2) = thetan(N);
        POWER(loop1,loop2) = v0*fextf;
        VARIATION1(loop1,loop2) = max(dthetan*L./thetan);
        VARIATION2(loop1,loop2) = max(dthetan*(1./lambda)./thetan);
        VARIATION3(loop1,loop2) = (max(thetan)-min(thetan))/mean(thetan);
        VARIATION4(loop1,loop2) = mean(dthetan*L./thetan);
        VARIATION5(loop1,loop2) = int_vn/(ksigman/(eta+etast))*(kad+dg)/(etast*Jactinf)/tauD;
        DIFFERENCE(loop1,loop2) = (Jactinf ...
            - dx/2*(sum(gamma(1:N-1).*thetan(1:N-1))+sum(gamma(2:N).*thetan(2:N))))/Jactinf*100;
        
        if loop1 == N1 && loop2 == N2 
            figure(21);
            subplot(3,1,1)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,gamma,x*1d-3,ksigman*thetan);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'\gamma (1/s)','fontsize',18) % left y-axis
            ylabel(hAx(2),'\sigma (Pa)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(3,1,2)
            [hAx,hLine1,hLine2] = plotyy(x*1d-3,thetan*1d3,x*1d-3,thetac*1d3);
            xlabel('x ({\mu}m)','fontsize',18)
            ylabel(hAx(1),'\theta_n ({\mu}M)','fontsize',18) % left y-axis
            ylabel(hAx(2),'\theta_c ({\mu}M)','fontsize',18) % right y-axis
            set(gca,'fontsize',18); set(hAx(2),'fontsize',18);
            set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
            vaxis = axis; vaxis(2) = L*1d-3; axis(vaxis);
            subplot(3,1,3)
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
    
    figure(103)
    plot(x*1d-3,thetan*1d3,'-','linewidth',2); hold on
    set(gca,'fontsize',18); 
    xlabel('x - x^b({\mu}m)','fontsize',18)
    ylabel('\theta_n ({\mu}M)','fontsize',18)
    xlim([0 L*1d-3])
    box off
    legend('\gamma = 5\times10^{-4} 1/s','\gamma = 1\times10^{-3} 1/s',...
        '\gamma = 5\times10^{-3} 1/s','\gamma = 1\times10^{-2} 1/s')

    figure(104)
    plot(x*1d-3,vn,'-','linewidth',2); hold on
    set(gca,'fontsize',18);
    xlabel('x - x^b ({\mu}m)','fontsize',18)
    ylabel('v_n (nm/s)','fontsize',18) 
    xlim([0 L*1d-3])
    box off
    legend('\gamma = 5\times10^{-4} 1/s','\gamma = 1\times10^{-3} 1/s',...
        '\gamma = 5\times10^{-3} 1/s','\gamma = 1\times10^{-2} 1/s')
    
%     figure(105)
%     plot(x*1d-3,dthetan*1d6,'-','linewidth',2); hold on
%     set(gca,'fontsize',18);
%     xlabel('x ({\mu}m)','fontsize',18)
%     ylabel('\theta_n{\prime} (mM/{\mu}m)','fontsize',18) 
%     xlim([0 L*1d-3])
%     box off
%     
%     figure(106)
%     sum2 = ksigman/(eta+etast)*thetan + Dtc*thetac;
%     plot(x*1d-3,sum2*1d-3,'-','linewidth',2); hold on
%     set(gca,'fontsize',18);
%     xlabel('x ({\mu}m)','fontsize',18)
%     ylabel('D_{\theta_n}\theta_n + D_{\theta_c}\theta_c ({\mu}m^2{\mu}M/s)','fontsize',18) 
%     xlim([0 L*1d-3])
%     box off
    
end


if Parameter1 > 0 && Parameter2 > 0
    [Xmesh,Ymesh] = meshgrid(Xm,Ym);
    
    figure(301)
    [~,hd] = contour(Xmesh,Ymesh,V0,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(V0)),max(max(V0)),25));
    title('v_0 (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    
    figure(302)
    subplot(2,2,1)
    [~,hd] = contour(Xmesh,Ymesh,V0,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(V0)),max(max(V0)),25));
    title('v_0 (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    subplot(2,2,2)
    [~,hd] = contour(Xmesh,Ymesh,VNF,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(VNF)),max(max(VNF)),25));
    title('v_n^f (nm/s)','fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    subplot(2,2,3)
    [~,hd] = contour(Xmesh,Ymesh,TNB,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(TNB)),max(max(TNB)),25));
    title('\theta_n^b ({\mu}M)','fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    subplot(2,2,4)
    [~,hd] = contour(Xmesh,Ymesh,TNF,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(TNF)),max(max(TNF)),25));
    title('\theta_n^f ({\mu}M)','fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    
    figure(303)
    [~,hd] = contour(Xmesh,Ymesh,DIFFERENCE,'fill','on'); colorbar;
    set(hd,'LevelList',linspace(min(min(DIFFERENCE)),max(max(DIFFERENCE)),25));
    title('(J_{actin}^f - \int_0^L{\gamma}\theta_ndx)/J_{actin}^f (%)',...
        'fontsize',18,'FontWeight','normal')
    AxesProperties20201003
    
    if Parameter2 == 9 
        figure(304)
        index = POWER < -inf;
        POWER(index) = NaN;
        [~,hd] = contour(Xmesh,Ymesh,POWER,'fill','on'); colorbar;
        set(hd,'LevelList',linspace(min(min(POWER)),max(max(POWER)),25));
        title('Power (nW/m^2)','fontsize',18,'FontWeight','normal')
        AxesProperties20201003
        
        
        TargetGamma = [5d-4, 8d-4,1d-3,2d-3];
        TargetGamma = [5d-3, 8d-3,1d-2,2d-2];
        nT = length(TargetGamma);
        for inT = 1:nT
            [~,idx] = min(abs(Ym-TargetGamma(inT)));
            figure(3043)
            vcell = V0(idx,:);
            power = POWER(idx,:);
            index = vcell >=0;
            plot(vcell(index),power(index),'-','linewidth',2); hold on
            set(gca,'fontsize',18);
            xlabel('v_0 (nm/s)')
            ylabel('Power (nW/m^2)')
            box off
%             legend('k_{\sigma_a} = 1\times10^{2} Pa/mM','k_{\sigma_a} = 5\times10^{2} Pa/mM',...
%                 'k_{\sigma_a} = 6\times10^{3} Pa/mM');
            figure(3044)
            vcell = V0(idx,:);
            power = POWER(idx,:);
            index = vcell >=0;
            plot(Xm(index),vcell(index),'-','linewidth',2); hold on
            set(gca,'fontsize',18);
            xlabel('f_{stall}^f (Pa)')
            ylabel('v_0 (nm/s)')
            box off


        end
    end
    
    if Parameter1 == 9
        
        figure(305)
        index = POWER < -inf;
        POWER(index) = NaN;
        [~,hd] = contour(Xmesh,Ymesh,POWER,'fill','on'); colorbar;
        set(hd,'LevelList',linspace(min(min(POWER)),max(max(POWER)),25));
        title('Power (nW/m^2)','fontsize',18,'FontWeight','normal')
        AxesProperties20201003
        
        Power = max(POWER);
        figure(306)
        plot(Xm,Power,'-','linewidth',2); hold on
        AxesProperties20201003
        ylabel('Max. Power (nW/m^2)')
        vaxis = axis;
        vaxis(3) = min(Power);
        vaxis(4) = max(Power);
        axis(vaxis);
        box off
    end
    
    
    
    
end

if Parameter1 == 0 && Parameter2 > 0
    figure(201)
    plot(Xm,V0,'-','linewidth',2); hold on
    AxesProperties20201003
    ylabel('v_0 (nm/s)')
    box off
    grid off
%     legend('k_{\sigma_n} = 1\times10^4 Pa/mM','k_{\sigma_n} = 2\times10^4 Pa/mM',...
%         'k_{\sigma_n} = 3\times10^4 Pa/mM','k_{\sigma_n} = 4\times10^4 Pa/mM')
    
    if Parameter2 == 9 
        figure(202)
        index = V0 >= 0;
        plot(V0(index),POWER(index),'-','linewidth',2); hold on
        set(gca,'fontsize',18);
        xlabel('v_0 (nm/s)')
        ylabel('Power (nW/m^2)')
        box off
    end
    
%     figure(2031)
%     plot(Xm,VARIATION1,'-','linewidth',2); hold on
%     AxesProperties20201003
%     ylabel('max(L\theta_n{\prime}/\theta_n)')
%     box off
% %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
% %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
%     
%     figure(2032)
%     plot(Xm,VARIATION2,'-','linewidth',2); hold on
%     AxesProperties20201003
%     ylabel('max(\theta_n{\prime}/(\lambda\theta_n))')
%     box off
% %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
% %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
%     
%     figure(2033)
%     plot(Xm,VARIATION3,'-','linewidth',2); hold on
%     AxesProperties20201003
%     ylabel('\Delta\theta_n/ mean(\theta_n)')
%     box off
% %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
% %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
%     
%     figure(2034)
%     plot(Xm,VARIATION4,'-','linewidth',2); hold on
%     AxesProperties20201003
%     ylabel('mean(L\theta_n{\prime}/\theta_n)')
%     box off
% %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
% %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
%     
%     figure(2035)
%     loglog(Xm,VARIATION5,'-','linewidth',2); hold on
%     AxesProperties20201003
%     box off
%     grid on
% %     legend('k_{\sigma_n} = 1\times10^2 (Pa/mM)','k_{\sigma_n} = 2\times10^2 (Pa/mM)',...
% %         'k_{\sigma_n} = 3\times10^2 (Pa/mM)','k_{\sigma_n} = 4\times10^2 (Pa/mM)')
    
end




