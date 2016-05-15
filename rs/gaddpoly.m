function z = gaddpoly (x, y)
% the less array index, the less power of element
% index=1 --> power=0
M = length(x);
N = length(y);
if M > N
    z = x;
else
    z = y;
end
for i=1:min(M,N)
    z(i) = bitxor(x(i), y(i));
end
