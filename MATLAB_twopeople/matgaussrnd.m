function [ X ] = matgaussrnd( M, U, V )
%MATGAUSSRND Matrix Gaussian random numbers

vecM = M(:);
S = kron(V,U);

Xvec = mvnrnd(vecM, S)';

X = reshape(Xvec, size(M));

end

