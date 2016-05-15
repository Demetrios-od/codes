function G = genConvPoly (g)
% g is a matrix of generators
% each row contains a generator in octal form
[M N] = size(g);
G = zeros(1,N*3);
pp = 1;
for col=1:N
    for bs=3:-1:1
        for row=1:M
            G(pp) = bitset(G(pp), M+1-row, bitget(g(row,col), bs));
        end
        pp = pp+1;
    end
end
while G(1) == 0
    G(1) = [];
end
