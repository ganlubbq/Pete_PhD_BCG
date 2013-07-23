clup

% model
K = 100;
ds = 2;
do = 2;
m = zeros(ds,1);
P = [1 -0.1; -0.1 1];
H = [1 0; 0 0.1];
R = [1 0.5; 0.5 1];

% data
w = mvnrnd(m',P)';
y = zeros(do,K);
for kk = 1:K
    y(:,kk) = mvnrnd((H*w)', R)';
end

% KF likelihood
kf_lhood = zeros(1,K);
kf_m = m;
kf_P = P;
for kk = 1:K
    kf_lhood(kk) = loggausspdf(y(:,kk), H*kf_m, H*kf_P*H'+R);
    [kf_m, kf_P] = kf_update(kf_m, kf_P, y(:,kk), H, R);
end

% Total likelihood
tmp = repmat({R}, K, 1);
tot_R = blkdiag(tmp{:});
tot_H = repmat(H, K, 1);
tot_lhood = loggausspdf(y(:), tot_H*m, tot_H*P*tot_H'+tot_R);