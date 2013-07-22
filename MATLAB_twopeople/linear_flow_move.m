function [ x, prob_ratio, drift] = linear_flow_move( t, t0, x0, m, P, y, H, R )
%linear_flow_move Calculate the movement produced by a linear flow for a
%set of matrixes. Optimal for linear Gaussian models. Analytic integral of
%linear_flow.

% Useful matrixes
HRH = H'*(R\H); HRH = (HRH+HRH')/2;

% OID approximations
Sigma0 = inv(inv(P)+t0*HRH);
Sigmat = inv(inv(P)+t*HRH);
mu0 = Sigma0*(t0*H'*(R\y)+P\m);
mut = Sigmat*(t*H'*(R\y)+P\m);
sqrtSigmat = sqrtm(Sigmat);

% Transformation
Gam = sqrtSigmat/sqrtm(Sigma0);
x = mut + Gam*(x0-mu0);

% Jacobian
prob_ratio = sqrt(det(Sigmat)/det(Sigma0));

% Final drift
drift = Sigmat*(H'/R)*( (y-H*mut)-0.5*H*(x-mut) );

end
