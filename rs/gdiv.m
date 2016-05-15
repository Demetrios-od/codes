function a = gdiv(x, y)
if y == 0
    a = NaN;
    return
end
if x == 0
    a = 0;
    return
end
lx = glog2(x);
ly = glog2(y);
la = lx-ly;
if la <= 0
    la = la+15;
end
a = gpow2(la);