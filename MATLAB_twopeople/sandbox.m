%% Correlation testing

C = [];
for ss = 1:4
    
    C(:,ss) = xcorr(observ(:,ss));
    
end

figure, plot(time, C(ceil(length(C)/2):end,:))
figure, hold on, plot(time, sum(C(ceil(length(C)/2):end,1:2),2)), plot(time, sum(C(ceil(length(C)/2):end,3:4),2), 'r')
figure, plot(time, sum(C(ceil(length(C)/2):end,:),2))

%%

