function a = gmul(x, y)
% works only within 4-digit Galois field (0,15)
if x==0 || y==0
    a = 0;
    return
end
lx = glog2(x);
ly = glog2(y);
la = mod(lx+ly, 15);
a = gpow2(la);