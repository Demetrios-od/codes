function a = gpow2(x)
if x >= 8
    x = mod(x,7);
end
if x == 0
    a = 1;
    return
end
base = [2 4 3 6 7 5 1];
a = base(x);