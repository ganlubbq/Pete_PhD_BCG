%% Correlation testing

K = 100;

C = [];
for ss = 1:model.num_sens
    
    C(:,ss) = xcorr(observ(ss,1:K));
    
end

figure, plot(time(1:K), C(ceil(length(C)/2):end,:))
figure, hold on, plot(time(1:K), sum(C(ceil(length(C)/2):end,1:2),2)), plot(time(1:K), sum(C(ceil(length(C)/2):end,3:4),2), 'r')
figure, plot(time(1:K), sum(C(ceil(length(C)/2):end,:),2))

%%

