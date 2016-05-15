function codes ()
clc
% parameters of Reed-Solomon (15,9) code
n = 15;    % words in the resulting block
k = 9;     % user-defined words to transmit
%m = 4;     % bits in data block - really is not need

% create a message for transmission (k positions)
%msg = round((2^m-1)*rand(1,k));   % random data
msg = [7 5 10 0 9 1 1 1 9];       % example from the book
%msg = [7 4 9 11 7 12 4 11 10];     % my example
disp('Original message:')
disp(msg)

%======================= ENCODING =======================
disp('ENCODING')
% create the generation polynom for RS code
%g = createGPolynom(n-k+1, m);
g = [12 10 12 3 9 7 1];   % for our initial data

% create the simplest generation matrix
G1 = zeros(k,n);
for i=1:k
    G1(i,:) = [zeros(1,i-1) g zeros(1,k-i)];
end
%G1

% create the generation matrix in systematic way
P = zeros(k, n-k);
for i=1:k
    P(i,:) = gdivpoly([zeros(1, n-i) 1], g);
end
G = [eye(k) P];

% create the parity-check matrix
H = [P' eye(n-k)];
% check: if G and H are correct, then obtain zero matrix
if ~any(gmulmatr(G, H'))
    s = 'OK';
else
    s = 'FAIL';
end
disp(['Parity-check matrix: ', s])

% calculate the code word using polynom g
% we'll build a systematic code, so now shift data to the right
msgext = [zeros(1,n-k) msg];

% divide extended msg to the generation polynom
% and obtain the encoded message
parity = gdivpoly(msgext, g);
cmsgp = [parity msg];

% check for errors
if ~any(gdivpoly(cmsgp, g))  % if empty matrix - then no errors
   s = 'OK';
else
    s = 'FAIL';
end
disp(['Coding by polynom: ', s])

% calculate the code word using non-systematic matrix G1
cmsgG1 = zeros(1,n);
for i=1:k
    cmsgG1 = gaddpoly(cmsgG1, gmulpoly(msg(i), G1(i,:)));
end
%cmsgG1
if ~any(gdivpoly(cmsgG1, g))   % should be empty matrix
    s = 'OK';
else
    s = 'FAIL';
end
disp(['Coding by non-systematic matrix: ', s])

% calculate the code word using systematic matrix G
cmsgG = zeros(1,n);
for i=1:k
    cmsgG = gaddpoly(cmsgG, gmulpoly(msg(i), G(i,:)));
end
%cmsgG

ee = gaddpoly(gdivpoly([zeros(1,n-k) cmsgG(1:k)], g), cmsgG(k+1:n)); % ???
if ~any(ee)
    s = 'OK';
else
    s = 'FAIL';
end
disp(['Coding by systematic matrix: ', s])

%======================= TRANSMITTING =======================
msgtr = cmsgp;   % choose message to transmit
% use a code word obtained by polynom

% add errors to the message
% ner = 2;                          % number of errors
% ev = round((2^m-1)*rand(1,ner));  % error values
% ep = 1+round((n-1)*rand(1,ner));  % error positions

% for testing - apply pre-defined errors
ev = [0 0 0 0 0 0 0 0 7 0 0 3 0 2 0];

% insert errors into message
%msgtr(ep) = ev;   % simply replace
msgtr = gaddpoly(msgtr, ev);

%======================= DECODING =======================
disp('DECODING')
ee = gaddpoly(gdivpoly([zeros(1,n-k) msgtr(n-k+1:n)], g), msgtr(1:n-k));
if ~any(ee)
    disp('There were no errors in the message.')
    disp('The message was:')
    disp(msgtr(n-k+1:n))
    return
end
% here we passed if errors were
disp('There were errors, start decoding')

% divide corrupted message to polynom
e = gdivpoly(msgtr, g);
S = zeros(size(e));
for i=1:length(e)
    S(i) = gcalcpoly(e, gpow2(i));
end
%S   % syndrome

% here we have to calculate locator polynom SIG
% use Massey algorithm (described in Morelos-Zaragoza's book, page 114)
sig = 1;
ro = [0 1];
le = 0;
for i=1:length(S)
    d = S(i);
    for j=1:le
        d = bitxor(d, gmul(sig(j+1), S(i-j)));
    end
    if d ~= 0
        signew = gaddpoly(sig, gmulpoly(d, ro));
        if le<i
            le = i-le;
            ro = zeros(1, length(sig));
            for j=1:length(sig)
                ro(j) = gdiv(sig(j), d);
            end
        end
        sig = signew;
    end
    ro = [0 ro];
end
sig
% I don't know why the given algorithm doesn't work
% The example from the book is passed well

sig = [1 6 5 4]   % these values must be for given data

W = gmulpoly(sig, S);
W(n-k:end) = [];

sig1 = sig(2:end);
for i=2:2:length(sig1)
    sig1(i) = 0;
end
%sig1   % derivative of SIG

v = zeros(1,length(sig));
for i=1:15    % find roots by checking each possible value
    v(i) = gcalcpoly(sig,i);
end
roots = find(~v);
pos = zeros(1, length(roots));
for i=1:length(roots)
    pos(i) = glog2(gdiv(1, roots(i)));
end
%pos   % error positions

evrx = zeros(1, length(roots));    % error values in received code word
for i=1:length(roots)
    evrx(i) = gdiv(gcalcpoly(W, roots(i)), gcalcpoly(sig1, roots(i)));
end
%evrx

msgtr(pos+1) = gaddpoly(msgtr(pos+1), evrx);
%msgtr

ee = gaddpoly(gdivpoly([zeros(1,n-k) msgtr(n-k+1:n)], g), msgtr(1:n-k));
if ~any(ee)
    disp('The message was:')
    disp(msgtr(n-k+1:n))
else
    disp('Decoding failed')
end
