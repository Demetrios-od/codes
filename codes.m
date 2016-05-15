function codes()
clc
u = [1 0 1 0 1 1 1 0 0 1 0 1];   % message: 12 bits
u = round(rand(1,12));
%u = [1 1 0 1 0 0];   % test
disp('Original message: ')
disp(u)

% matrix of generators in octal form
g = [1 3 3;
     1 7 1;
     1 6 5];   % my code
%g = [7; 5];    % test code from book
[M N] = size(g);   % M - number of generators

%========================== ENCODING ==========================

% make generation polynom - an element of the generation matrix
% here we need generators represented as a matrix with octal numbers
% G = zeros(1,N*3);
% pp = 1;
% for col=1:N
%     for bs=3:-1:1
%         for row=1:M
%             G(pp) = bitset(G(pp), M+1-row, bitget(g(row,col), bs));
%         end
%         pp = pp+1;
%     end
% end
% while G(1) == 0
%     G(1) = [];
% end
% G

% but it'll be better if generators would be represented
% as binaries (decimals)

% convert generators from octal to binary (decimal) form
g1 = zeros(1,M);
for row=1:M
    for col=1:N
        g1(row) = g1(row)+g(row,col)*8^(N-col);
    end
end

% define how many bits maximum in each generator
v = 1;   % how many bits
f = true;
while f
    v = bitshift(v,1);
    f = false;
    for i=1:M
        f = f || g1(i)>v;
    end
end

% make generation polynom
G = 0;
i = 1;
while v>1
    G(i) = 0;
    for j=1:M
        G(i) = bitset(G(i), M+1-j, bitget(g1(j), i));
    end
    i = i+1;
    v = bitshift(v,-1);
end
G = fliplr(G);   % = [7 3 7 6 1 4 7];

emsg = gmulpoly(u,G);     % encoding message

% convert encoded message into array with 0 and 1
msgtr = zeros(1, length(emsg)*M);
for i=0:length(emsg)-1
    for j=1:M
        msgtr(i*M+j) = bitget(emsg(i+1), M+1-j);
    end
end
% disp('Message without errors')
% disp(msgtr)

%====================== INSERTING ERRORS ======================

m = length(G)-1;  % number of digits in coder's register

err = zeros(1, length(msgtr));
ner = 2;    % number of errors
ep = 1+round((length(msgtr)-1)*rand(1,ner));  % error positions
err(ep) = 1;
disp(['Number of errors: ', int2str(length(find(err)))])
msgtr = gaddpoly(msgtr, err);

% convert into BPSK symbols with noise
msgtr = round(4*rand(1,length(msgtr))) - 4*msgtr;

% disp('Message with errors')
% disp(msgtr)

%========================== DECODING ==========================

NN = 2^m;         % number of nodes (states)

% calculate trellis matrix
T = -ones(NN);    % empty value = -1
for row=0:NN-1    % transferring from this state...
    for bin=0:1   % input bit can be 0 or 1, only two cases
        col = bitset(bitshift(row, -1), m, bin);   % ...to this state
        s = bitset(row, m+1, bin);   % state of coder's register
        out = 0;
        for i=1:M    % for each generator
            b = 0;
            for j=1:m+1   % calculating each bit after current generator
                b = bitxor(b, bitand(bitget(s,j), bitget(g1(i),j)));
            end
            out = bitset(out, M+1-i, b);
        end
        T(row+1, col+1) = out;   % element of the trellis matrix
    end
end
%T

NM = -ones(1,NN);      % node metrics
NM(1) = 0;
for i=1:NN
    Npath{i} = [];   % saved path ended in each node
end
for i=1:M:length(msgtr)   % for each group of bits (symbols) at the input
    for j=1:NN
        NMvar{j,1} = [];     % variants for node metric
        NMvar{j,2} = [];
    end
    for j=0:NN-1       % for each node
        if NM(j+1) ~= -1   % skip non-initialized nodes
            conns = T(j+1,:);   % row of connections from the trellis matrix
            for ic=1:NN         % for each possible connection
                if conns(ic) ~= -1    % skip unused connections
                    % here: ic - output node; conns(ic) - output symbols at
                    % branch
                    metric = calcMetric(msgtr(i:i+M-1), conns(ic), true);
                    NMvar{ic,1} = [NMvar{ic,1}, metric];
                    NMvar{ic,2} = [NMvar{ic,2}, j];   % which node we're coming from
                end
            end
        end
    end
    % now all nodes are considered
    Npathnew = Npath;
    NMnew = NM;
    for j=1:NN      % for each node
        v = NMvar{j,1};     % metrics of current branches
        pn = NMvar{j,2};    % previous nodes
        if ~isempty(v)
            % add metrics of previous nodes and current branches,
            % and get minimun; define position of min(max) element
            [NMnew(j), vp] = max(v + NM(pn+1));   % max for soft decision
            Npathnew{j} = [Npath{pn(vp)+1}, bitget(j-1,m)];
        end
    end
    Npath = Npathnew;
    NM = NMnew;
end
out = Npath{1};   % because our coder always returns to zero state
out(length(u)+1:end) = [];
disp('Decoded message:')
disp(out)
end


function m = calcMetric(ip, br, isSoftDecision)
% ip (input) - array from the input sequence (1 and 0 for hard decision)
% br (branch) - symbol at the trellis branch (a number with some digits)
M = length(ip);
m = 0;

if isSoftDecision
    for i=1:M
        m = m + ip(i)*(~bitget(br, M+1-i));
    end
else    % hard decision
    for i=1:M
        m = m + bitxor(ip(i), bitget(br, M+1-i));
    end
end
end
