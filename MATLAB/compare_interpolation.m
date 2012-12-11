clup

load template_beat.mat

set_model;
set_algo;

model.dw = 30;
template = template';

t = (0:0.001:2)';
t_in = (0:model.dw-1)/model.fs;

% Sinc
H = heartbeat_interpolation(algo, model, t, 0);
s_sinc = H*template;

% Spline
s_spline = interp1(t_in, template, t, 'spline', 0);

% Cubic
s_cubic = interp1(t_in, template, t, 'cubic', 0);

% Grid
n = 0:model.dw-1;
grid = (t*model.fs*ones(1,model.dw)-ones(length(t),1)*n);
lt1 = abs(grid)<1;
lt2 = abs(grid)<2;
lt3 = abs(grid)<3;

% Square-cosine
interp_mat = zeros(size(grid));
interp_mat(lt1) = cos(pi*grid(lt1)/2).^2;
s_cosine = interp_mat*template;

% Windowed-sinc
interp_mat = zeros(size(grid));
wind = zeros(size(grid));
wind(lt3) = 1;
interp_mat(lt3) = sinc(grid(lt3)).*wind(lt3);
s_windsinc = interp_mat*template;

% % sinc-polynomial
% lto = abs(grid)<1;
% interp_mat = zeros(size(grid));
% interp_mat(lto) = 1 - (pi*grid(lto)).^2/6 + (pi*grid(lto)).^4/120 - (pi*grid(lto)).^6/factorial(7);
% s_sincpoly = interp_mat*template;

figure, hold on
plot(t, s_sinc, 'g')
plot(t, s_spline, 'b')
plot(t, s_cubic, 'c')
plot(t, s_cosine, 'r')
plot(t, s_windsinc, 'm')
% plot(t, s_sincpoly, 'm')


%%

clear all

load template_beat.mat

set_algo;
set_model;

model.dw = 30;
template = template';

tau = (0:0.001:5)';
t = (0:100)'/model.fs;
y = [template' template' template' zeros(1,11)];

loglhood = zeros(size(tau));
dloglhood = zeros(size(tau));
% ccor = zeros(size(tau));
% reg = zeros(size(tau));
for kk = 1:length(tau)
    H = heartbeat_interpolation(algo, model, t, tau(kk));
    
    n = 0:model.dw-1;
    grid = (t*model.fs*ones(1,model.dw)-ones(length(t),1)*n) - tau(kk)*model.fs;
    H1 = -sincd(grid,1)*model.fs;
    
    loglhood(kk) = -0.5*sum((y' - H*template).^2)/model.y_obs_vr;
    dloglhood(kk) = (y*H1*template/model.y_obs_vr);
%     ccor(kk) = -2*sum(y*H*template);
%     reg(kk) = template'*H'*H*template;
end

figure, hold on, plot(t,y)
figure, hold on, plot(tau, loglhood, 'k')
figure, hold on, plot(tau(2:end), diff(loglhood)/0.001, 'r'), plot(tau, dloglhood, 'b')

% test laplace approximation
tau0 = 0.97;
n = 0:model.dw-1;
grid = (t*model.fs*ones(1,model.dw)-ones(length(t),1)*n) - tau0*model.fs;
H1 = -sincd(grid,1)*model.fs;
H2 = sincd(grid,2)*(model.fs^2);

vr = -1./(y*H2*template/model.y_obs_vr);
mn = tau0 + vr*(y*H1*template/model.y_obs_vr);
figure, plot(tau, (-(tau-mn).^2)./vr, 'r')