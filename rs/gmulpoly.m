function z = gmulpoly(x, y)
% the less array index, the less power of element
% index=1 --> power=0
M = length(x)-1;
N = length(y)-1;
z = zeros(1, M+N);
for i=0:M+N
    z(i+1) = 0;
    for j=max(0,i-N):min(i,M)
        z(i+1) = bitxor(z(i+1), gmul(x(j+1), y(i-j+1)));
    end
end
