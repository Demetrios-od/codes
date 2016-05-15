function C = gmulmatr(A, B)
sa = size(A);
sb = size(B);
if sa(2) ~= sb(1)
   error('MULMATR error: Matrices must match'); 
end
C = zeros(sa(1), sb(2));
for i=1:sa(1)
    for j=1:sb(2)
        for k=1:sa(2)
            C(i,j) = bitxor(C(i,j), gmul(A(i,k), B(k,j)));
        end
    end
end