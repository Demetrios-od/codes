clc
% g = [1 3 3;
%      1 7 1;
%      1 6 5];   % my code
% %g = [7; 5];    % test code from book
% [M N] = size(g)
% 
% g1 = zeros(1,M);
% for row=1:M
%     for col=1:N
%         g1(row) = g1(row)+g(row,col)*8^(N-col);
%     end
% end
% g1
% 
% v = 1;
% f = true;
% while f
%     v = bitshift(v,1);
%     f = false;
%     for i=1:M
%         f = f || g1(i)>v;
%     end
% end
% v

ip = [0 1];
br = 2;
M = length(ip);
m = 0;
for i=1:M
    m = m + bitxor(ip(i), bitget(br, M+1-i));
end
m