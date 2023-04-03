% This is part of the orginal codes used in the following paper:
% http://www.molbiolcell.org/cgi/doi/10.1091/mbc.E22-10-0494
% On the role of myosin-induced actin depolymerization during cell migration
% If you have questions, feel free to contact Dr. Yizeng Li.

% Analytical solution for linear Jactinf in thetac

clear
clc


FadMode = 1;
% FadMode = 1: Fad = kad*v0
% FadMode = 2: Fad = kad*thetac(1)*v0

JactinMode = 1;
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
Theta  = 0.3;           % (mM) reference total G- and F-actin, \int_0^L (thetan + thetac)dx/L

Mc = 5d-3;              % (mM) reference value of free myosin
Mn = 5d-3;              % (mM) reference value of

if JactinMode == 1
    Jactinf0 = 30;          % (nm/s) Jactin^f = Jactinf0*thetac^f
elseif JactinMode == 2
    Jactinf0 = 6;           % (nm mM/s) Jactinf = Jactinf0*thetac^f/(thetacc + thetac^f)
    thetacc  = 0.2d-3;         % (mM) Critical value for actin polymerization
end

if GammaMode == 1
    gamma0 = 5d-4;          % (1/s) constant rate of actin depolymerization
    % gamma0 = 1d-3;
    % gamma0 = 5d-3;
    % gamma0 = 1d-2;
elseif GammaMode == 2
    gamma0 = 1d-4;      % (1/s) gamma = gamma0 + gammmaa*sigma_a
    gammaa = 2d-3;      % (1/Pa/s)
elseif GammaMode == 3
    gamma0 = 5d-4;      % (1/s) gamma = gamma0 + gammaa*sigma_a/(sigma_ac + sigma_a)
    gammaa = 5d-3;      % (1/s)
    sigma_ac = 1;       % (Pa)
end

ka = 1d-0;     % (1/s mM) ka*m_c*theta_n
ko = 1.5d-0;   % (1/s) ko*m_n 

ksigman = 100*1d2;          % (Pa /mM) Coefficient of passive actin network pressure
ksigmaa = 4d2;          % (Pa /mM) Coefficient of active actin network contraction
 
etast = 100*1d-4;           % (Pa s/nm^2/mM)
eta   = 1d-8;           % (Pa s/nm^2/mM)
dg    = 1d-6;           % (Pa s/nm) coefficient of hydraulic resistance
if FadMode == 1
    kad = 100*3d-1;           % (Pa s/nm) adhesive force, Fad^b = kad*v0
elseif FadMode == 2
    kad = 1d0;          % (Pa s/mM/nm) adhesive force, Fad^b = kad*thetan^b*v0
end

Dtc = 1.d7;             % (nm^2/s) diffusion constant for theta_c
Dmc = 1.d6;             % (nm^2/s) diffusion constant for m_c
Dmn = 1.d5;             % (nm^2/s) diffusion constant for m_n

N = 101;
gamma = logspace(-4,0,N); % (1/s)


lambda = sqrt((etast)/ksigman)*gamma.^(1/2);

v0 = Jactinf0*Theta*ksigman*(exp(lambda*L)+exp(-lambda*L)-2)...
    ./((kad+dg)*((ksigman*lambda/etast + Jactinf0./lambda/L.*(1-ksigman/Dtc/etast))...
    .*(exp(lambda*L)-exp(-lambda*L)) + Jactinf0*ksigman/Dtc/etast*(exp(lambda*L)+exp(-lambda*L)))...
    +2*Jactinf0*Theta*etast./lambda.*(1-exp(-lambda*L))...
    +Jactinf0*Theta*eta./lambda.*(exp(lambda*L)-exp(-lambda*L)));


delta2 = 2*etast*v0./ksigman./lambda;


figure(1)
semilogx(gamma,v0,'-','linewidth',2); hold on
set(gca,'fontsize',15);
xlabel('\gamma (1/s)')
ylabel('v_0 (nm/s)')
box off
legend('k_{\sigma_n} = 1\times10^4 Pa/mM','k_{\sigma_n} = 2\times10^4 Pa/mM',...
    'k_{\sigma_n} = 3\times10^4 Pa/mM','k_{\sigma_n} = 4\times10^4 Pa/mM')

figure(2)
semilogx(gamma,delta2,'-','linewidth',2); hold on
set(gca,'fontsize',15);
xlabel('\gamma (1/s)')
ylabel('\delta_2 (nm/s)')
box off
legend('k_{\sigma_n} = 1\times10^4 Pa/mM','k_{\sigma_n} = 2\times10^4 Pa/mM',...
    'k_{\sigma_n} = 3\times10^4 Pa/mM','k_{\sigma_n} = 4\times10^4 Pa/mM')
