%% Correlation testing

K = 100;

C = [];
for ss = 1:model.num_sens
    
    C(:,ss) = xcorr(observ(ss,1:K));
    
end

figure, plot(time(1:K), C(ceil(length(C)/2):end,:))
figure, hold on, plot(time(1:K), sum(C(ceil(length(C)/2):end,1:2),2)), plot(time(1:K), sum(C(ceil(length(C)/2):end,3:4),2), 'r')
figure, plot(time(1:K), sum(C(ceil(length(C)/2):end,:),2))

%% Likelihood optimisation testing

dt = 0.001;
tau_rng = -2:dt:3;

func=zeros(1,length(tau_rng));
grad=zeros(1,length(tau_rng));
hess=zeros(1,length(tau_rng));

for ll=1:length(tau_rng), [func(ll), grad(ll), hess(ll)] = log_lhood(algo, model, time, observ, beat, p_idx, b_idx, wf_mn, wf_vr, tau_rng(ll)); end

figure
plot(tau_rng, func)

figure, hold on
plot(tau_rng, grad)
plot(tau_rng, [0 func(3:end)-func(1:end-2) 0]/(2*dt), 'r')
% plot(tau_rng, [0 0 -(1/12)*func(5:end)+(2/3)*func(4:end-1)-(2/3)*func(2:end-3)+(1/12)*func(1:end-4) 0 0]/dt, 'r')
% plot(tau_rng, [0 diff(func)]/dt, 'g')