function a = gdivpoly(x, y)
% divides polynoms: x/y
% polynom y must have the largest (last) coefficient = 1
% x(1), y(1) - the lowest powers (=0)
% result - module after division
N = length(y);
if length(x)<N
    a = x;
    return
end
x = fliplr(x);
y = fliplr(y);
while length(x) >= N
    r = x(1);
    for i=1:N
        x(i) = bitxor(x(i), gmul(y(i), r));
    end
    while ~isempty(x) && x(1) == 0
        x = x(2:end);
    end
end
a = fliplr(x);
