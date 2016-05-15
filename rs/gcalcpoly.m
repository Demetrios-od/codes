function a = gcalcpoly(p, x)
% p - polynom to calculate
% x - value to substitute
a = 0;
lx = glog2(x);
for i=0:length(p)-1
    p2 = gpow2(mod(lx*i,15));
    elem = gmul(p(i+1), p2);
    a = bitxor(a, elem);
end
