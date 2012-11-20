function [ pdf ] = loggampdf( x,a,b )
%LOGGAMPDF Finds the log-density of a point from a gamma distribution

pdf = - a*log(b) - gammaln(a) + (a-1)*log(x) - x/b;

end

