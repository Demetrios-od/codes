function a = glog2(x)
while x >= 8
    x = mod(x,7);
end
if x == 0
    a = NaN;
    return
end
base = [7 1 3 2 6 4 5];
a = base(x);
