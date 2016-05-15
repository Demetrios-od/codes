function a = createGPolynom(D)
% D - the max power of the polynom
if D == 1
    a = 1;
    return
end
for i=1:D-1
    a = gmulpoly([gpow2(i) 1], createGPolynom(D-1));
end
