% This is part of the orginal codes used in the following paper:
% http://www.molbiolcell.org/cgi/doi/10.1091/mbc.E22-10-0494
% On the role of myosin-induced actin depolymerization during cell migration
% If you have questions, feel free to contact Dr. Yizeng Li.

if Parameter1 == 1
    set(gca,'fontsize',18,'yscale','log');
    if JactinMode == 1
        ylabel('J_a (nm/s)','fontsize',18)
    elseif JactinMode == 2
        ylabel('J_a (nm mM/s)','fontsize',18)
    end
elseif Parameter1 == 2
    set(gca,'fontsize',18,'yscale','log');
    ylabel('\eta_{st} (Pa s/{\mu}m^2/mM)','fontsize',18)
elseif Parameter1 == 3
    set(gca,'fontsize',18,'yscale','log');
    ylabel('\gamma (1/s)','fontsize',18)
elseif Parameter1 == 4
    set(gca,'fontsize',18,'yscale','log');
    if FadMode == 1
        ylabel('k_{ad} (Pa s/{\mu}m)','fontsize',18)
    elseif FadMode == 2
        ylabel('k_{ad} (Pa s//mM/{\mu}m)','fontsize',18)
    end
elseif Parameter1 == 5
    set(gca,'fontsize',18,'yscale','log');
    ylabel('d_g (Pa s/{\mu}m)','fontsize',18)
elseif Parameter1 == 6
    set(gca,'fontsize',18,'yscale','log');
    ylabel('\eta (Pa s/{\mu}m^2/mM)','fontsize',18);
elseif Parameter1 == 7
    set(gca,'fontsize',18,'yscale','log');
    ylabel('D_{\theta_c} ({\mu}m^2/s)','fontsize',18);
elseif Parameter1 == 8
    set(gca,'fontsize',18,'yscale','log');
    ylabel('k_{\sigma_n} (Pa/mM)','fontsize',18)
elseif Parameter1 == 9
    set(gca,'fontsize',18);
    ylabel('f_{ext}^f (Pa)','fontsize',18)
end

if Parameter2 == 1
    set(gca,'fontsize',18,'xscale','log');
    if JactinMode == 1
        xlabel('J_a (nm/s)','fontsize',18)
    elseif JactinMode == 2
        xlabel('J_a (nm mM/s)','fontsize',18)
    end
elseif Parameter2 == 2
    set(gca,'fontsize',18,'xscale','log');
    xlabel('\eta_{st} (Pa s/{\mu}m^2/mM)','fontsize',18)
elseif Parameter2 == 3
    set(gca,'fontsize',18,'xscale','log');
    xlabel('\gamma (1/s)','fontsize',18)
elseif Parameter2 == 4
    set(gca,'fontsize',18,'xscale','log');
    if FadMode == 1
        xlabel('k_{ad} (Pa s/{\mu}m)','fontsize',18)
    elseif FadMode == 2
        xlabel('k_{ad} (Pa s//mM/{\mu}m)','fontsize',18)
    end
elseif Parameter2 == 5
    set(gca,'fontsize',18,'xscale','log');
    xlabel('d_g (Pa s/{\mu}m)','fontsize',18)
elseif Parameter2 == 6
    set(gca,'fontsize',18,'xscale','log');
    xlabel('\eta (Pa s/{\mu}m^2/mM)','fontsize',18);
elseif Parameter2 == 7
    set(gca,'fontsize',18,'xscale','log');
    xlabel('D_{\theta_c} ({\mu}m^2/s)','fontsize',18);
elseif Parameter2 == 8
    set(gca,'fontsize',18,'xscale','log');
    xlabel('k_{\sigma_n} (Pa/mM)','fontsize',18)
elseif Parameter2 == 9
    set(gca,'fontsize',18);
    xlabel('f_{ext}^f (Pa)','fontsize',18)
end


if Parameter1 > 0 && Parameter2 > 0   
    vaxis = axis;
    vaxis(1) = min(min(Xmesh));
    vaxis(2) = max(max(Xmesh));
    vaxis(3) = min(min(Ymesh));
    vaxis(4) = max(max(Ymesh));
    axis(vaxis);
    
elseif Parameter1 == 0 && Parameter2 > 0  
    
end