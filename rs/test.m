clc
clear
% disp('testing glog2')
% n = 20;
% z = zeros(2,N);
% for i=0:N
%     z(:,i+1) = [i, glog2(i)];
% end
% z

% disp('-------------------------------------------------------')
% disp('testing gpow2')
% n = 6;
% z = zeros(2,n+1);
% for i=0:n
%     z(:,i+1) = [i, gpow2(i)];
% end
% z

% disp('-------------------------------------------------------')
% disp('testing gmul')
% for x=0:15
%     for y=0:15
%         if gmul(x,y) ~= gmulold(x,y,4)
%             disp('Error occurs in GMUL')
%             x
%             y
%         end
%     end
% end
% N = 15;
% s = zeros(2,15-N);
% for i=N:15
%     s(:,i-N+1) = [i, gmul(i,N)];
% end
% s

% disp('-------------------------------------------------------')
% disp('testing gdiv')
% for x=0:15
%     for y=1:15
%         z = gmul(x,y);
%         if gdiv(z,y) ~= x
%             disp('Error occurs in GDIV')
%             x
%             y
%         end
%     end
% end

% disp('-------------------------------------------------------')
% disp('testing gaddpoly')
% x = [2 3 8 12]
% y = [9 14 3]
% gaddpoly(x,y)

% disp('-------------------------------------------------------')
% disp('testing gmulpoly')
% x = [1 12 0 4]
% y = [9 5]
% gmulpoly(x,y)

% disp('-------------------------------------------------------')
% disp('testing gdivpoly')
% x = [1 12 0 4 7 15 6 2]
% y = [9  10 1]   % last position (largest power) must be 1
% gdivpoly(x,y)

disp('-------------------------------------------------------')
disp('testing gcalcpoly')
% p = [12 10 12 3 9 7 1]
% s = zeros(1,6);
% for i=1:6
%     s(i) = gcalcpoly(p,gpow2(i));
% end
% s

p = [1 6 5 4];
a = zeros(2,length(p));
for i=1:15
    a(:,i) = [i, gcalcpoly(p,i)];
end
a
x = find(~a(2,:))
for i=1:length(x)
    x(i) = glog2(gdiv(1, x(i)));   % THESE ARE ERROR LOCATIONS !!!
end
x


% disp('-------------------------------------------------------')
% disp('testing createGPolynom')
% createGPolynom(7)

% disp('-------------------------------------------------------')
% disp('testing gmulmatr')
% x = [ 3  7 12;
%      14  1  0;
%       5  6  2]
% y = [ 8  1 10;
%       3  9 11;
%      15 12  4]
% gmulmatr(x,y)

disp('-------------------------------------------------------')
disp('test is finished')
