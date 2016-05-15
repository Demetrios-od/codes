function a = gpow2(x)
if x >= 16
    x = mod(x,15);
end
if x == 0
    a = 1;
    return
end
base = [2 4 8 3 6 12 11 5 10 7 14 15 13 9 1];
a = base(x);