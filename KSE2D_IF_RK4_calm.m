function [output] = KSE2D_IF_RK4_calm(N,T,lambda, epsilon, type, iniData)
close all;

%% (0) Default Parameters ================================================

if ~(exist('N','var')) % Number of x gridpoints 
    N = 128;
end
if ~(exist('T','var')) % Final time
    T = 20;
end
if ~(exist('lambda','var')) % Instability parameter
    lambda = 4.1;
end
if ~(exist('epsilon','var')) % Calming parameter
    epsilon = 1e-4;
end
if ~(exist('type','var')) % Calming type (set = 1, 2, or 3)
    type = 2;
end
if ~(exist('iniData','var')) % Initial data type (set = 1 or 2)
    iniData = 1;
end

% Parameter settings
p.N         = N;
p.T         = T;
p.C         = 0.5;
p.type      = type;
p.eps       = epsilon;
p.lambda    = lambda;

%% (1) Parameters/Settings/Update times ===================================
setting.norms_on        = 0; % 1 = compute norms over time 
setting.spec_on         = 0; % 1 = compute spectrum when figure displays
setting.save_on         = 0; % 1 = Saves figures to folder whenever display is updated
setting.gif_on          = 0; % 1 = Saves gif of program
setting.display_on      = 0; % 1 = Displays solution and spectrum for eye-norm (See (8D))
setting.convergence_on  = 1; % 1 = Compute error between KSE solution and calmed KSE solution

movie_speed = 0.2; % delay time between frames

setting.update_time = p.T/100; % Update figures every h seconds.
fig_dt = setting.update_time;

Markers=['o','+','*','s','d','v','^','<','>','p','x']; % To use for marking lines in plots

seed = 0;
rng(seed);

%% (2) Name tags ==========================================================
tag = ['_lam_' num2str(p.lambda) '_N_' num2str(p.N) '_type_' num2str(p.type) '_eps_' num2str(p.eps)];
base_name = ['cKSE', tag]; base_name = strrep(base_name,'.','p');
spec_name = ['spec', tag]; spec_name = strrep(spec_name,'.','p');

movname = ['gifs\' base_name '.gif'];
calm_tag = ['cKSE', tag];

if (setting.gif_on == 1)
    movname = [base_name '.gif'];
end

if (setting.save_on == 1)
    save_dir = makeFolder(mfilename, base_name, 'saved_figures');
end

%% (3) Initial Data =======================================================
N = p.N;
L_x = 2*pi; dx = L_x/N; x_0 = -pi; x = x_0 + dx*(0:(N-1));
L_y = 2*pi; dy = L_y/N; y_0 = -pi; y = y_0 + dy*(0:(N-1));

u_0 = zeros(N,N);
v_0 = zeros(N,N);


for n = 1: N
    for m = 1: N
        if iniData == 1
            u_0(m,n) = (cos(x(m) + y(n)) + cos(x(m)));
            v_0(m,n) = (cos(x(m) + y(n)) + cos(y(n)));
        elseif iniData == 2
            u_0(m,n) = 4*(cos(x(m) + y(n)) + sin(3*x(m)));
            v_0(m,n) = 4*(cos(x(m) + y(n)) + cos(4*y(n)));
        end
    end
end

u_KSE_hat  = fftn(u_0);
v_KSE_hat  = fftn(v_0);

u_calm_hat = fftn(u_0);
v_calm_hat = fftn(v_0);


%==========================================================================
%==========================================================================

%% (4) Time mesh ==========================================================
t0 = 0; T = p.T;

CFL_adv = 0.5;

max_adv_vel = p.C*1.7*p.lambda^2*max(max(abs([u_0, v_0])));

dt = (dx/max_adv_vel)*CFL_adv;

t = t0:dt:T;

%==========================================================================
%==========================================================================

%% (5) Derivatives in Fourier Space =======================================
p.kx = [0:(N/2), ((-N/2)+1):-1]*(2*pi/L_x);
p.ikx = 1i*p.kx;
p.ikx((N/2)+1) = 0; % Kill Nyquist frequency in x for odd-order derivatives
p.ky = [0:(N/2), ((-N/2)+1):-1]*(2*pi/L_y);
p.iky = 1i*p.ky;
p.iky((N/2)+1) = 0; % Kill Nyquist frequency in y for odd-order derivatives
p.ksq = zeros(N,N);
p.dealias_modes = [1+ ceil(N/3): N+1 - ceil(N/3)];
for j = 1:N % Iterating rows
    for k = 1:N % Iterating columns
        p.ksq(j,k) = p.kx(j)^2 + p.ky(k)^2;
    end
end
%==========================================================================
%==========================================================================

%% (6) Spectral Operators =================================================

p.L_full = exp(dt*(p.lambda*p.ksq - (p.ksq).^2)  );
p.L_half = exp(0.5*dt*(p.lambda*p.ksq - (p.ksq).^2));

%==========================================================================
%==========================================================================

%% (7) Spectrum Plotting ==================================================
% Display the spectrum of the solution over time.

spectrum_KSE  = zeros(1, floor(sqrt(2)*(N) ));
spectrum_calm = zeros(1, floor(sqrt(2)*(N) ));

if setting.spec_on == 1
    %     if setting.convergence_on == 1
    %         spectrum_error = zeros(1, floor(sqrt(2)*(N) ));
    %     end
    for m = 1: (1+N/2) % Iterate over the unique values of |k|.
        for n = 1: (1+N/2)
            wave_mode = floor(sqrt(m^2 + n^2));
            spectrum_KSE(wave_mode) = spectrum_KSE(wave_mode) + abs(u_KSE_hat(m,n))^2 + abs(v_KSE_hat(m,n))^2;
            spectrum_calm(wave_mode) = spectrum_calm(wave_mode) + abs(u_calm_hat(m,n))^2 + abs(v_calm_hat(m,n))^2;
        end
    end
    spectrum_KSE = (sqrt(spectrum_KSE))/N^2;
    spectrum_calm = (sqrt(spectrum_calm))/N^2;
    
    spec_check = spectrum_KSE(floor(1+(N/3))); % Measure highest value of spectrum near Nyquist frequency
    
end

%% (8) Figure Display =====================================================

fig = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(sprintf('lambda = %g, time = 0', p.lambda));

KSE = subplot(1,3,1);
p_KSE = pcolor(x, y, ((ifftn(u_KSE_hat, 'symmetric').^2 + ifftn(v_KSE_hat, 'symmetric').^2).^(0.5)));
colorbar;
set(p_KSE, 'Edgecolor', 'None');
shading interp; lighting gouraud; axis('tight'); axis('square');
title('KSE');
% drawnow;

cKSE = subplot(1,3,2);
p_cKSE = pcolor(x, y, ((ifftn(u_calm_hat, 'symmetric').^2 + ifftn(v_calm_hat, 'symmetric').^2).^(0.5)));
colorbar;
set(p_cKSE, 'Edgecolor', 'None');
shading interp; lighting gouraud; axis('tight'); axis('square');
title(sprintf('cKSE type %g, epsilon = %f',p.type, p.eps ));
% drawnow;

spec = subplot(1,3,3);

mi = 1;
loglog( (1:(1+N/2)), spectrum_KSE((1:(1+N/2))), 'k', 'DisplayName', 'KSE'); hold on;
loglog( (1:(1+N/2)), spectrum_calm((1:(1+N/2))), 'b-x', 'DisplayName', 'cKSE');
xline(1+(N/3), '--', 'linewidth', 2);
yline(eps, ':', 'linewidth', 2);
legend('KSE', 'cKSE', 'N/3', '2.2204e-16','location', 'southoutside', 'fontsize', 14);

axis('tight');
axis('square');
xlabel('Wave number k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',10);
hold off;
drawnow;

% Start gif
if setting.gif_on == 1
    movie_frame = getframe(fig);
    theIm = frame2im(movie_frame);
    [imind,cm] = rgb2ind(theIm,256);
    imwrite(imind, cm, movname, 'gif', 'DelayTime', movie_speed, 'Loopcount', inf);
end

%% (9) Norm setup =========================================================

p.parseval = sqrt(L_x*L_y)/(N^2);
p.parseval_sq = p.parseval^2;

if setting.convergence_on == 1
    output.LinfL2_error      = 0;
    output.L2H2_error_sq     = 0;
    output.LinfLinf_error    = 0;
end

%% (10) Time loop =========================================================

for j = 2:(length(t))
    % Note: Pre-RK4 step, we have u_hat ~ u_hat(t(j-1))
    
    %% (A) Runge-Kutta Scheme =============================================
    
    if setting.convergence_on == 1
        % Convergence of cKSE to KSE
        % Stage 1:
        % Update KSE and calm KSE solutions at the same time,
        % compute error in norms.
        
        [norms, K1, H1, k1, h1] = ...
            rhs_error(p, u_KSE_hat, v_KSE_hat, u_calm_hat, v_calm_hat);
        
        % Update remaining stages of KSE solutions
        %         [u_new, v_new, n] = rhs(norms_on, p, u_hat, v_hat)
        [K2, H2, ~] = rhs_KSE(0, p,   p.L_half.*(u_KSE_hat + 0.5*dt.*K1), ...
            (                     p.L_half.*(v_KSE_hat + 0.5*dt.*H1)));
        [K3, H3, ~] = rhs_KSE(0, p,   p.L_half.*u_KSE_hat + 0.5*dt*K2, ...
            (                     p.L_half.*v_KSE_hat + 0.5*dt*H2) );
        [K4, H4, ~] = rhs_KSE(0, p,   p.L_full.*u_KSE_hat + dt*p.L_half.*K3,...
            (                     p.L_full.*v_KSE_hat + dt*p.L_half.*H3));
        % Update remaining stages of calm KSE solutions
        [k2, h2, ~] = rhs_cKSE(0, p,   p.L_half.*(u_calm_hat + 0.5*dt.*k1), ...
            (                     p.L_half.*(v_calm_hat + 0.5*dt.*h1)));
        [k3, h3, ~] = rhs_cKSE(0, p,   p.L_half.*u_calm_hat + 0.5*dt*k2, ...
            (                     p.L_half.*v_calm_hat + 0.5*dt*h2) );
        [k4, h4, ~] = rhs_cKSE(0, p,   p.L_full.*u_calm_hat + dt*p.L_half.*k3,...
            (                     p.L_full.*v_calm_hat + dt*p.L_half.*h3));
        % Apply final RK4 step to both solutions
        u_calm_hat     = p.L_full.*u_calm_hat + (dt/6)*( p.L_full.*k1 + 2*p.L_half.*(k2 + k3) + k4);
        v_calm_hat     = p.L_full.*v_calm_hat + (dt/6)*( p.L_full.*h1 + 2*p.L_half.*(h2 + h3) + h4);
        
        u_KSE_hat = p.L_full.*u_KSE_hat + (dt/6)*( p.L_full.*K1 + 2*p.L_half.*(K2 + K3) + K4);
        v_KSE_hat = p.L_full.*v_KSE_hat + (dt/6)*( p.L_full.*H1 + 2*p.L_half.*(H2 + H3) + H4);
        % =========================================================================================
    else % Else, just run cKSE
        % =========================================================================================
        
        [k1, h1, output] = rhs_cKSE(setting.norms_on, p, u_calm_hat, v_calm_hat);
        
        [k2, h2, ~]      = rhs_cKSE(0, p,   p.L_half.*(u_calm_hat + 0.5*dt.*k1), ...
            (                     p.L_half.*(v_calm_hat + 0.5*dt.*h1)));
        [k3, h3, ~]      = rhs_cKSE(0, p,   p.L_half.*u_calm_hat + 0.5*dt*k2, ...
            (                     p.L_half.*v_calm_hat + 0.5*dt*h2) );
        [k4, h4, ~]      = rhs_cKSE(0, p,   p.L_full.*u_calm_hat + dt*p.L_half.*k3,...
            (                     p.L_full.*v_calm_hat + dt*p.L_half.*h3));
        
        u_calm_hat = p.L_full.*u_calm_hat + (dt/6)*( p.L_full.*k1 + 2*p.L_half.*(k2 + k3) + k4);
        v_calm_hat = p.L_full.*v_calm_hat + (dt/6)*( p.L_full.*h1 + 2*p.L_half.*(h2 + h3) + h4);
        
    end
    
    % Note: It is after the RK4 step that u_hat ~ u_hat(t(j))
    %======================================================================
    
    %% (B) Update norms ===================================================
    
    if setting.norms_on == 1
        % Compute Linf
        Linf = max( max(abs([ ifftn(u_calm_hat, 'symmetric'), ifftn(v_calm_hat, 'symmetric') ]) ) );
        % Compute LinfLinf
        Linf_MAX = max( [Linf_MAX, Linf] );
    end
    
    if setting.convergence_on == 1
        output.LinfL2_error      = max( output.LinfL2_error, sqrt( norms.L2_error_sq  ) );
        
        Linf_error = max( norm( ifftn(u_calm_hat - u_KSE_hat , 'symmetric'),'inf'), ...
            norm( ifftn(v_calm_hat - v_KSE_hat , 'symmetric'),'inf'));
        output.LinfLinf_error    = max( output.LinfLinf_error, Linf_error);
        
        output.L2H2_error_sq     = output.L2H2_error_sq + (norms.L2_error_sq + norms.L2_lap_error_sq)*dt;
    end
    %======================================================================
    
    %% (C) Display/Save figures ====================================================
    if setting.display_on == 1
        if (fig_dt - t(j) < 0.5*dt )  % Display figures
            
            sgtitle(sprintf('lambda = %g, time = %1.4f', p.lambda, t(j) ));
            
            % Compute spectrum
            if (setting.spec_on == 1)
                spectrum_KSE(:) = 0; spectrum_calm(:) = 0;
                for m = 1: N % Iterate over the unique values of |k|.
                    for n = 1: N
                        wave_mode = floor(sqrt(m^2 + n^2));
                        spectrum_KSE(wave_mode) = spectrum_KSE(wave_mode) + abs(u_KSE_hat(m,n))^2 + abs(v_KSE_hat(m,n))^2;
                        spectrum_calm(wave_mode) = spectrum_calm(wave_mode) + abs(u_calm_hat(m,n))^2 + abs(v_calm_hat(m,n))^2;
                    end
                end
                spectrum_KSE = (sqrt(spectrum_KSE))/N^2;
                spectrum_calm = (sqrt(spectrum_calm))/N^2;
                spec_check = max([spectrum_KSE(floor(1+(N/3))), spectrum_calm(floor(1+(N/3)))]);
                
                if spec_check > 1e-8
                    warning('Spectrum at Nyquist frequency is above single-precision.')
                end
                
                spec;
                cla(spec);
                
                loglog( (1:(1+N/2)), spectrum_KSE((1:(1+N/2))),  'k', 'DisplayName', 'KSE'); hold on;
                loglog( (1:(1+N/2)), spectrum_calm((1:(1+N/2))), 'b-x', 'DisplayName', 'cKSE');
                xline(1+(N/3), '--', 'linewidth', 2);
                yline(eps, ':', 'linewidth', 2);
                legend('KSE', 'cKSE', 'N/3', '2.2204e-16','location', 'southoutside', 'fontsize', 14);
                axis('tight');
                axis('square');
            end
            
            KSE;   set(p_KSE,   'cdata', ((ifftn( u_KSE_hat, 'symmetric').^2 + ifftn(v_KSE_hat, 'symmetric').^2).^(0.5)));
            cKSE;  set(p_cKSE,  'cdata', ((ifftn(u_calm_hat, 'symmetric').^2 + ifftn(v_calm_hat, 'symmetric').^2).^(0.5)));
            drawnow;
            
            % Reset position of figure before saving pictures.
            fig.OuterPosition = [0 0 1 1];
            
            % Make gif
            if (setting.gif_on == 1)
                movie_frame = getframe(fig);
                theIm = frame2im(movie_frame);
                [imind,cm] = rgb2ind(theIm, 256);
                imwrite(imind, cm, movname, 'gif', 'DelayTime', movie_speed, 'WriteMode','append');
            end
            
            % Save picture
            if (setting.save_on == 1)
            saveFigure(save_dir, base_name, t(j));
            end
                        
            fig_dt = fig_dt + setting.update_time;
        end
    end
    
end % End Time loop =======================================================

%% (11) Update terminal data

if (setting.spec_on == 1)
    % Compute spectrum
    spectrum_KSE(:) = 0; spectrum_calm(:) = 0;
    for m = 1: N % Iterate over the unique values of |k|.
        for n = 1: N
            wave_mode = floor(sqrt(m^2 + n^2));
            spectrum_KSE(wave_mode) = spectrum_KSE(wave_mode) + abs(u_KSE_hat(m,n))^2 + abs(v_KSE_hat(m,n))^2;
            spectrum_calm(wave_mode) = spectrum_calm(wave_mode) + abs(u_calm_hat(m,n))^2 + abs(v_calm_hat(m,n))^2;
        end
    end
    spectrum_KSE = (sqrt(spectrum_KSE))/N^2;
    spectrum_calm = (sqrt(spectrum_calm))/N^2;
    spec_check = max([spectrum_KSE(floor(1+(N/3))), spectrum_calm(floor(1+(N/3)))]);
    
    
    if setting.display_on == 1
%         spec;
        cla(spec);
        
        loglog( (1:(1+N/2)), spectrum_KSE((1:(1+N/2))),  'k', 'DisplayName', 'KSE'); hold on;
        loglog( (1:(1+N/2)), spectrum_calm((1:(1+N/2))), 'b-x', 'DisplayName', 'cKSE');
        xline(1+(N/3), '--', 'linewidth', 2);
        yline(eps, ':', 'linewidth', 2);
        legend('KSE', 'cKSE', 'N/3', '2.2204e-16','location', 'eastoutside', 'fontsize', 14);
    end
end

if setting.display_on == 1
    sgtitle(sprintf('lambda = %g, time = %1.4f', p.lambda, t(j) ));
    
    KSE;   
    set(p_KSE,   'cdata', ((ifftn(u_KSE_hat, 'symmetric').^2 + ifftn(v_calm_hat, 'symmetric').^2).^(0.5)));
    cKSE;  
    set(p_cKSE,  'cdata', ((ifftn(u_calm_hat, 'symmetric').^2 + ifftn(v_calm_hat, 'symmetric').^2).^(0.5)));
    drawnow;
end

% Reset position of figure before saving pictures.
fig.OuterPosition = [0 0 1 1];

% Make gif
if (setting.gif_on == 1)
    movie_frame = getframe(fig);
    theIm = frame2im(movie_frame);
    [imind,cm] = rgb2ind(theIm, 256);
    imwrite(imind, cm, movname, 'gif', 'DelayTime', movie_speed, 'WriteMode','append');
end

% Take picture
if (setting.save_on == 1)
        saveFigure(save_dir, base_name, T);
end

end % End function KSE2D_IF_RK4_calm


%% (I) Right hand side error computation =====================================
function [norms, u_KSE_new, v_KSE_new, u_calm_new, v_calm_new] = ...
    rhs_error(p, u_KSE_hat, v_KSE_hat, u_calm_hat, v_calm_hat)
%   % N(u) = -eta_eps(u)u_x - eta_eps(v)u_y
%   % N(v) = -eta_eps(u)v_x - eta_eps(v)v_y

%% (A) Dealias modes =========================================================
u_calm_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for u in x
u_calm_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for u in y
v_calm_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for v in x
v_calm_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for v in y

u_KSE_hat(p.dealias_modes, :)  = 0; % Eliminates higher frequencies for u in x
u_KSE_hat(:, p.dealias_modes)  = 0; % Eliminates higher frequencies for u in y
v_KSE_hat(p.dealias_modes, :)  = 0; % Eliminates higher frequencies for v in x
v_KSE_hat(:, p.dealias_modes)  = 0; % Eliminates higher frequencies for v in y

%% (B) Compute derivatives ===================================================
% For calm solution
u_calm_x_hat    = bsxfun(@times, p.ikx.', u_calm_hat);
u_calm_y_hat    = bsxfun(@times, p.iky, u_calm_hat);
v_calm_x_hat    = bsxfun(@times, p.ikx.', v_calm_hat);
v_calm_y_hat    = bsxfun(@times, p.iky, v_calm_hat);
u_calm_x        = ifftn(u_calm_x_hat, 'symmetric');
u_calm_y        = ifftn(u_calm_y_hat, 'symmetric');
v_calm_x        = ifftn(v_calm_x_hat, 'symmetric');
v_calm_y        = ifftn(v_calm_y_hat, 'symmetric');
% For KSE solution
u_KSE_x_hat     = bsxfun(@times, p.ikx.', u_KSE_hat);
u_KSE_y_hat     = bsxfun(@times, p.iky, u_KSE_hat);
v_KSE_x_hat     = bsxfun(@times, p.ikx.', v_KSE_hat);
v_KSE_y_hat     = bsxfun(@times, p.iky, v_KSE_hat);
u_KSE           = ifftn(u_KSE_hat,   'symmetric');
u_KSE_x         = ifftn(u_KSE_x_hat, 'symmetric');
u_KSE_y         = ifftn(u_KSE_y_hat, 'symmetric');
v_KSE           = ifftn(v_KSE_hat,   'symmetric');
v_KSE_x         = ifftn(v_KSE_x_hat, 'symmetric');
v_KSE_y         = ifftn(v_KSE_y_hat, 'symmetric');

%% (C) Compute norms =========================================================

% Computes || u_hat ||^2_(L2) and || LAP(u_hat) ||^2_(L2)
norms.L2_error_sq     = p.parseval_sq*(norm(u_KSE_hat - u_calm_hat, 'fro')^2 + ...
    (                              norm(v_KSE_hat - v_calm_hat, 'fro')^2));
norms.L2_lap_error_sq = p.parseval_sq*(norm(-p.ksq.*(u_KSE_hat - u_calm_hat), 'fro')^2 + ...
    (                              norm(-p.ksq.*(v_KSE_hat - v_calm_hat), 'fro')^2));

%% (D) Initialize calming function ===========================================
if (p.type == 0)        % Type 0: Just KSE
    eta_eps_u = ifftn(u_calm_hat, 'symmetric');
    eta_eps_v = ifftn(v_calm_hat, 'symmetric');
elseif (p.type == 1)    % cKSE Type 1
    eta_eps_u = ifftn(u_calm_hat, 'symmetric')./( 1 + p.eps*abs(ifftn(u_calm_hat, 'symmetric')));
    eta_eps_v = ifftn(v_calm_hat, 'symmetric')./( 1 + p.eps*abs(ifftn(v_calm_hat, 'symmetric')));
elseif (p.type == 2)    % cKSE Type 2
    eta_eps_u = ifftn(u_calm_hat, 'symmetric')./( 1 + (p.eps*abs(ifftn(u_calm_hat, 'symmetric'))).^2 );
    eta_eps_v = ifftn(v_calm_hat, 'symmetric')./( 1 + (p.eps*abs(ifftn(v_calm_hat, 'symmetric'))).^2 );
elseif (p.type == 3)    % cKSE Type 3
    eta_eps_u = (1/p.eps)*atan(p.eps*ifftn(u_calm_hat, 'symmetric'));
    eta_eps_v = (1/p.eps)*atan(p.eps*ifftn(v_calm_hat, 'symmetric'));
end

%% (E) Update solutions ======================================================
u_calm_new = -fftn(  eta_eps_u.*u_calm_x + ...
    (                eta_eps_v.*u_calm_y));
v_calm_new = -fftn(  eta_eps_u.*v_calm_x + ...
    (                eta_eps_v.*v_calm_y));
u_KSE_new  = -fftn(  u_KSE.*u_KSE_x + ...
    (                v_KSE.*u_KSE_y));
v_KSE_new  = -fftn(  u_KSE.*v_KSE_x + ...
    (                v_KSE.*v_KSE_y));
end % End rhs_error

%% (II) Calm KSE Right hand side =============================================
function [u_new, v_new, n] = rhs_cKSE(norms_on, p, u_hat, v_hat)
%   % N(u) = -eta_eps(u)u_x - eta_eps(v)u_y
%   % N(v) = -eta_eps(u)v_x - eta_eps(v)v_y

u_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for u in x
u_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for u in y
v_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for v in x
v_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for v in y

if (p.type == 0)        % Type 0: Just KSE
    eta_eps_u = ifftn(u_hat, 'symmetric');
    eta_eps_v = ifftn(v_hat, 'symmetric');
elseif (p.type == 1)    % cKSE Type 1
    eta_eps_u = ifftn(u_hat, 'symmetric')./( 1 + p.eps*abs(ifftn(u_hat, 'symmetric')));
    eta_eps_v = ifftn(v_hat, 'symmetric')./( 1 + p.eps*abs(ifftn(v_hat, 'symmetric')));
elseif (p.type == 2)    % cKSE Type 2
    eta_eps_u = ifftn(u_hat, 'symmetric')./( 1 + (p.eps*abs(ifftn(u_hat, 'symmetric'))).^2 );
    eta_eps_v = ifftn(v_hat, 'symmetric')./( 1 + (p.eps*abs(ifftn(v_hat, 'symmetric'))).^2 );
elseif (p.type == 3)    % cKSE Type 3
    eta_eps_u = (1/p.eps)*atan(p.eps*ifftn(u_hat, 'symmetric'));
    eta_eps_v = (1/p.eps)*atan(p.eps*ifftn(v_hat, 'symmetric'));
end

% Take derivatives
u_x_hat = bsxfun(@times, p.ikx.', u_hat);
u_y_hat = bsxfun(@times, p.iky, u_hat);
v_x_hat = bsxfun(@times, p.ikx.', v_hat);
v_y_hat = bsxfun(@times, p.iky, v_hat);

u_x = ifftn(u_x_hat, 'symmetric');
u_y = ifftn(u_y_hat, 'symmetric');
v_x = ifftn(v_x_hat, 'symmetric');
v_y = ifftn(v_y_hat, 'symmetric');

% Update
u_new = -fftn(  eta_eps_u.*u_x + ...
    (           eta_eps_v.*u_y));

v_new = -fftn(  eta_eps_u.*v_x + ...
    (           eta_eps_v.*v_y));

% Compute norms of the gradient
if (norms_on == 1)
    n.L2_grad = sqrt(norm(u_x_hat,'fro')^2 + norm(u_y_hat,'fro')^2 + norm(v_x_hat,'fro')^2 + norm(v_y_hat,'fro')^2)*p.parseval;
    n.L_inf_grad = max(max(abs([u_x, u_y, v_x, v_y]) ));
else
    n.L2_grad = 0;
    n.L_inf_grad = 0;
end

end % End rhs
%==========================================================================
%==========================================================================

%% (III) KSE Right hand side =================================================
function [u_new, v_new, n] = rhs_KSE(norms_on, p, u_hat, v_hat)
%   % N(u) = -eta_eps(u)u_x - eta_eps(v)u_y
%   % N(v) = -eta_eps(u)v_x - eta_eps(v)v_y

u_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for u in x
u_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for u in y
v_hat(p.dealias_modes, :) = 0; % Eliminates higher frequencies for v in x
v_hat(:, p.dealias_modes) = 0; % Eliminates higher frequencies for v in y

u = ifftn(u_hat, 'symmetric');
v = ifftn(v_hat, 'symmetric');


% Take derivatives
u_x_hat = bsxfun(@times, p.ikx.', u_hat);
u_y_hat = bsxfun(@times, p.iky, u_hat);
v_x_hat = bsxfun(@times, p.ikx.', v_hat);
v_y_hat = bsxfun(@times, p.iky, v_hat);

u_x = ifftn(u_x_hat, 'symmetric');
u_y = ifftn(u_y_hat, 'symmetric');
v_x = ifftn(v_x_hat, 'symmetric');
v_y = ifftn(v_y_hat, 'symmetric');

% Update
u_new = -fftn(  u.*u_x + ...
    (           v.*u_y));

v_new = -fftn(  u.*v_x + ...
    (           v.*v_y));


% Compute norms of the gradient
if (norms_on == 1)
    n.L2_grad = sqrt(norm(u_x_hat,'fro')^2 + norm(u_y_hat,'fro')^2 + norm(v_x_hat,'fro')^2 + norm(v_y_hat,'fro')^2)*p.parseval;
    n.L_inf_grad = max(max(abs([u_x, u_y, v_x, v_y]) ));
else
    n.L2_grad = 0;
    n.L_inf_grad = 0;
end

end % End rhs