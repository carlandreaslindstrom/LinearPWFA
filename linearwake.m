function [max_g] = linearwake(sigx_um, sigz_um, Qb_pC, p_mbar, visu, r_res, z_res)

    % define colormap
    function map = cmap(mn,mx)
        % determine relative position of the zero
        if nargin < 1; zero = 0.5; else; zero = -mn/(mx-mn); end;
        % COLORS
        if zero == 0; map = interp1([0, 0.5, 1],[255 255 255; 70 172 255; 0 51 139]/256,linspace(0, 1, 500));
        elseif zero == 1; map = interp1([0, 0.5, 1],[255,148,0; 255 197 124; 255 255 255]/256,linspace(0, 1, 500));
        else; map = interp1([0 mean([0,zero]) zero mean([zero,1]) 1],[255,148,0; 255 197 124; 255 255 255; 70 172 255; 0 51 139]/256,linspace(0, 1, 500)); 
        end
    end


    % DEFINE PARAMETERS
    % visualize by default
    if ~exist('visu','var'); visu = true; end;
    
    % "fundamental" constants
    SI_kb = 1.38064852e-23; % [J/K] Boltzmann constant
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_eps0 = 8.854187817e-12; % [F/m] permittivity of free space
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_c = 299792458; % [m/s] speed of light
    
    % plasma density
    p0 = p_mbar * 100; % [Pa] gas pressure
    T = 300; % [K] room temperature
    Nelectron = 1; % number of elecrons ionized
    n0 = Nelectron * p0 / (SI_kb * T); % [m^-3] 
    
    % beam density
    sigx = sigx_um * 1e-6; % [m] transverse size
    sigz = sigz_um * 1e-6; % [m] rms bunch length
    Qb = Qb_pC * 1e-12; % [C] bunch charge
    nb = -(Qb/SI_e)/(sqrt(2*pi)^3 * sigx^2 * sigz);
    
    % plasma wavenumber
    k_p = sqrt(n0*SI_e^2/(SI_eps0*SI_me*SI_c^2)); % [1/m]
    lambda_p = 2*pi/k_p; % [m]
    
    % box size and auto-resolution
    r_edge = 5*sigx;
    z_edge = 5*sigz;
    if ~exist('r_res','var'); r_res = 51; end;
    if ~exist('z_res','var')
        res_per_lambdap = 10;
        min_z_res = 51;
        z_res = max(res_per_lambdap*(z_edge*2)/lambda_p, min_z_res);
    end;
    r_res = floor(r_res/2)*2+1; % make odd
    z_res = floor(z_res/2)*2+1; % make odd
    
    % linear density perturbation (nb < n0)
    rs = linspace(-r_edge, r_edge, r_res);
    zs = linspace(-z_edge, z_edge, z_res);
    nbs = -(Qb/SI_e) * normpdf(rs, 0, sigx)'*normpdf(zs, 0, sigz)/sqrt(2*pi*sigx^2);
    
    
    % CALCULATE DENSITY PERTURBATION
    if visu; disp('Calculating densities...'); end
    dns = zeros(numel(rs), numel(zs));
    for i = 1:numel(rs)
        for j = numel(zs):(-1):1
            z = zs(j); % lower bound
            zprimes = zs(zs >= z);
            argument = nbs(i,j:numel(zs)) .* sin(k_p*(z-zprimes));
            dns(i,j) =  - k_p * trapz(zprimes, argument', 1);
        end
        
        % display progress bar
        if visu && mod(i,50)==0; fprintf([num2str(round(i/numel(rs)*100)) '%% ']); end
    end
    if visu; disp('... Done!'); end
    
    
    % CALCULATE FIELDS
    if visu; disp('Calculating fields...'); end
    zs_f = zs(2:(end-1));
    rs_f = rs(2:(end-1));
    Ezs = zeros(numel(rs_f), numel(zs_f));
    Ers = zeros(numel(rs_f), numel(zs_f));
    ddndz = (dns(2:(end-1),3:end)-dns(2:(end-1),1:(end-2)))/(2*mean(diff(zs)));
    ddndr = (dns(3:end,2:(end-1))-dns(1:(end-2),2:(end-1)))/(2*mean(diff(rs)));
    
    rmask = (rs_f >= -mean(diff(rs_f))/2);
    for i = 1:numel(rs_f) 
        r = rs_f(i);
        rs_prime = rs_f(rmask);
        r_gt = max(rs_prime, abs(r)); % r_>
        r_lt = min(rs_prime, abs(r)); % r_<
        zArg = rs_prime .* besselk(0, k_p*r_gt) .* besseli(0, k_p*r_lt);
        rArg = rs_prime .* besselk(1, k_p*r_gt) .* besseli(1, k_p*r_lt);
        zArg(isnan(zArg)) = 0; % fix for 0 argument in bessel function
        rArg(isnan(rArg)) = 0; % fix for 0 argument in bessel function
        for j = 1:numel(zs_f)
            Ezs(i,j) = SI_e/SI_eps0 * trapz(rs_prime, ddndz(rmask,j)' .* zArg);
            Ers(i,j) = - sign(r) * SI_e/SI_eps0 * trapz(rs_prime, ddndr(rmask,j)' .* rArg);
        end
        
        % display progress bar
        if visu && mod(i,50)==0; fprintf([num2str(round(i/numel(rs_f)*100)) '%% ']); end
    end
    if visu; disp('... Done!'); end
    
    
    % on axis field gradient
    jmid = round(numel(rs_f)/2);
    gs = 1/SI_c*(Ers(jmid+1,:) - Ers(jmid-1,:))/(2*mean(diff(rs_f)));
    max_g = max(abs(gs)); % find maximum gradient
    
    
    % PLOTS
    if visu
    
        figure(1);
        set(gcf,'color','w');
        colormap(cmap());

        % beam density
        subplot(4,2,1);
        imagesc(zs*1e6, rs*1e6, nbs/1e6);
        xlabel('z [\mum]');
        ylabel('r [\mum]');
        title(['Beam density (n_b = ' num2str(abs(nb/1e6),'%2.2g') '/cc)']);
        c = colorbar;
        ylabel(c,'Beam density [1/cc]');
        caxis([min(min(nbs)), -min(min(nbs))]/1e6);

        % plasma density
        subplot(4,2,2);
        imagesc(zs*1e6, rs*1e6, dns/1e6);
        xlabel('z [\mum]');
        ylabel('r [\mum]');
        title(['Plasma density (n_0 = ' num2str(abs(n0/1e6),'%2.2g') '/cc)']);
        c = colorbar;
        ylabel(c,'Density perturbation [1/cc]');
        caxis([min(min(dns)), -min(min(dns))]/1e6);

        % longitudinal field
        subplot(4,2,3);
        imagesc(zs_f*1e6, rs_f*1e6, Ezs*1e-9);
        xlabel('z [\mum]');
        ylabel('r [\mum]');
        title('Longitudinal field');
        c = colorbar;
        ylabel(c,'Accelerating gradient [GeV/m]');

        % transverse field
        subplot(4,2,4);
        imagesc(zs_f*1e6, rs_f*1e6, Ers*1e-9);
        xlabel('z [\mum]');
        ylabel('r [\mum]');
        title('Transverse field');
        c = colorbar;
        ylabel(c,'Transverse force [GeV/m]');

        % equivalent B-field
        subplot(4,2,5:6);
        [~, mid_beam_ind] = min(abs(zs_f));
        [~, one_sig_ind] = min(abs(zs_f - sigz_um*1e-6));
        plot(rs_f*1e6, Ers(:,mid_beam_ind)/SI_c, ...
            rs_f*1e6, Ers(:,one_sig_ind)/SI_c, 'LineWidth',2);
        xlim([min(rs_f), max(rs_f)]*1e6);
        xlabel('r [\mum]');
        ylabel('B_{equiv} [T]');
        title('B-field (equivalent) vs. radius');
        legend('Mid-beam', 'At +1\sigma_z');

        % magnetic field gradient vs. z
        subplot(4,2,7:8);
        yyaxis left;
        plot(zs_f*1e6, gs,'LineWidth',2); 
        xlabel('z [\mum]');
        ylabel('g [T/m]');
        title(['On-axis magnetic field gradient vs. z (max ' num2str(max_g,'%g') ' T/m)']);
        yyaxis right;
        plot(zs_f*1e6, Qb/SI_e*normpdf(zs_f, 0, sigz)/(2*pi*sigx^2),'LineWidth',2); 
        ylabel('dQ/dz [C/m]');
        
    end
    
end
