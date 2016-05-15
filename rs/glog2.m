function a = glog2(x)
while x >= 16
    x = mod(x,15);
end
if x == 0
    a = NaN;
    return
end
base = [15 1 4 2 8 5 10 3 14 9 7 6 13 11 12];
a = base(x);
