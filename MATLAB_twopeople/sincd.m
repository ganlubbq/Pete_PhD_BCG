function [ y ] = sincd( t, d )
%SINCD first two derivatives of sinc(x)

if nargin<2
    d=1;
else
    assert((d==1)||(d==2))
end

x = pi*t;
z = x==0;
x(z) = 1;

y = zeros(size(t));

if d==1
    y(~z) = (x(~z).*cos(x(~z))-sin(x(~z)))./(x(~z).^2);
    y(z) = 0;
    
elseif d==2
    y(~z) = ((2-x(~z).^2).*sin(x(~z))-2*x(~z).*cos(x(~z)))./(x(~z).^3);
    y(z) = -1/3;
    
end

end

