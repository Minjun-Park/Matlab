% Minjun Park, CAAM 210, 11/01/15
%
% bridgedriver.m
%
% solve and display a loaded bridge using Gaussian Elimination
%
% usage: bridgedriver(Ea,W,nos)
%
% where: Ea = fiber modulus-times-area
% W = vehicle weight
% nos = number of stages in bridge
%
% example: bridgedriver(1,.01,8);
%
function bridgedriver
close all %close everything before it begins.
Ea=1; %Set thickness (product of Youngs modulus and cross sectional
W=0.05; %Weight of the car
nos=8; %number of stages
bridge(Ea,W,nos) %call bridge
end
function bridge(Ea,W,nos)
close all
m = 5*nos - 5; % number of fibers
n = nos*2; % number of nodes
dof = 2*n; % number of degrees of freedom
s = 1/sqrt(2);
A = zeros(m,dof);
A(1,1) = 1; % build 1st corner
A(2,[3 4]) = [s s];
for i=1:nos-2 % build remaining stages
 hl = 4*(i-1) + 1; % index hor-dof of bottom left node
 A(3+(i-1)*5,[hl+1 hl+3]) = [-1 1];
 A(4+(i-1)*5,[hl+2 hl+3 hl+4 hl+5]) = [-s s s -s];
 A(5+(i-1)*5,[hl+2 hl+6]) = [-1 1];
 A(6+(i-1)*5,[hl hl+1 hl+6 hl+7]) = [-s -s s s];
 A(7+(i-1)*5,[hl hl+4]) = [-1 1];
end
hl = hl + 4;
A(m-2,[hl+1 hl+3]) = [-1 1];
A(m-1,[hl+2 hl+3]) = [-s s];
A(m,hl) = -1;
for i = 0: nos-2
A(16+3*i, 17+3*i)=1; % this is first stylin bridge
A(23+7*i, 16+7*i)=1; % this is first styline bridge
end
for i = 1:nos-2
    startfiber = 5*nos-4+i; % start of the stylin fiber
    A(startfiber, [dof-3 dof-2 3+4*(i-1) 4+4*(i-1)]) = [cos(atan(2/i)) -sin(atan(2/i)) -cos(atan(2/i)) sin(atan(2/i)) ];
    % this forms the rest of tall fibers
end

figure(1) % visualize the adjacencies
spy(A)
xlabel('degree of freedom','fontsize',16)
ylabel('fiber','fontsize',16)
title(['Nonzero adjacencies, nos = ' num2str(nos)],'fontsize',16)
xc = zeros(7*n-5,2);
yc = zeros(7*n-5,2);
xc([1 2],:) = [0 1; 0 1]; % left end
yc([1 2],:) = [0 0; 0 1];
for i=0:nos-3 % construct the coordinates of the undeformed net
 xc(3+i*5:7+i*5,:) = [i+1 i+1; i+1 i+2; i+1 i+2; i+1 i+2; i+1 i+2];

 yc(3+i*5:7+i*5,:) = [0 1; 1 0; 1 1; 0 1; 0 0];
end
xc([m-2:m],:) = [nos-1 nos-1; nos-1 nos; nos-1 nos];
yc([m-2:m],:) = [0 1; 1 0; 0 0];
xc([m+1],:) = [0 0]; % tallest fiber
yc([m+1],:) = [0 3]; % tallest fiber
xc(end,:)= [nos nos]; %last fiber(x coordinate)
yc(end,:)= [0 3]; %last fiber(y coordinate)
for i = 1: nos-1
   xc(m+1+i,:) = [i 0]; %left fibers of pylons(x-coordinate)
   yc(m+1+i,:) = [1 3]; %left fibers of pylons(y-coordinate)
end
for i = 1:nos-1
    xc(end-i,:) = [nos-i nos];%right fibers of pylons(x-coordinate)
    yc(end-i,:) = [1 3];%right fibers of pylons(y-coordinate)
end


figure(2)
line(xc',yc','linewidth',2) % draw the undeformed net
hold on
L = ones(m,1); % fiber lengths
L(2) = 1/s;
for j=1:nos-2,
 L(4+(j-1)*5) = 1/s;
 L(6+(j-1)*5) = 1/s;
end
L(m-1) = 1/s;
k = Ea./L; % the fiber stiffnesses
K = diag(k);
S = A'*K*A; % the mother stiffness matrix
f = zeros(dof,1); % the load on the bridge
for j = 1:nos-1,
 f(2+(j-1)*4) = -W; % downward force at subsequent nodes
 x = gauss(S,f); % using Gaussian elimination
 %x = S\f; % solve for the displacements % replace with gaussian
 X = x(1:2:dof);
 Y = x(2:2:dof);
 dx(1:2,:) = xc(1:2,:) + [0 X(1); 0 X(2)];
 dy(1:2,:) = yc(1:2,:) + [0 Y(1); 0 Y(2)];
 
 
% for i=0:nos-3 % construct the coordinates of the undeformed net
% xc(3+i*5:7+i*5,:) = [i+1 i+1; i+1 i+2; i+1 i+2; i+1 i+2; i+1 i+2];
%
% yc(3+i*5:7+i*5,:) = [0 1; 1 0; 1 1; 0 1; 0 0];
%end
 for i=0:nos-3 % get the deformed coordinates
 dx(3+i*5:7+i*5,:) = xc(3+i*5:7+i*5,:) + [X(1+2*i) X(2+2*i)
 X(2+2*i) X(3+2*i)
 X(2+2*i) X(4+2*i)
 X(1+2*i) X(4+2*i)
 X(1+2*i) X(3+2*i)];
 dy(3+5*i:7+5*i,:) = yc(3+i*5:7+i*5,:) + [Y(1+2*i) Y(2+2*i);
 Y(2+2*i) Y(3+2*i)
 Y(2+2*i) Y(4+2*i)
 Y(1+2*i) Y(4+2*i)
 Y(1+2*i) Y(3+2*i)];
 end
 dx(m+1,:) = xc(m+1,:) + [0 X(m+1)]; %deformed for pylon
 dy(m+1,:) = yc(m+1,:) + [0 Y(m+1)]; %deformed for pylon
 dx(end,:) = xc(end,:) + [0 X(end)]; %last fiber deformed(x coordinate)
 dy(end,:)=  dy(end,:) + [0 Y(end)]; %last fiber deformed(y coordinate)
 for i = 1: nos-1
 dx(m+1+i,:) = xc(m+1+i,:) + X(i); %left fibers of pylons deformed(x-coordinate)
 dy(m+1+i,:) = yc(m+1+i,:) + [1 3]; %left fibers of pylons deformed(y-coordinate)
 end
 for i = 1:nos-1
 dx(end-i,:) = xc(end-i,:) + [nos-i nos];%right fibers of pylons deformed(x-coordinate)
 dy(end-i,:) = yc(end-i,:) + [1 3];%right fibers of pylons deformed(y-coordinate)
 end
 
 dx(m-2:m,:) = xc([m-2:m],:) + [X(end-1) X(end)
 X(end) 0
 X(end-1) 0];
 dy(m-2:m,:) = yc([m-2:m],:) + [Y(end-1) Y(end)
 Y(end) 0
 Y(end-1) 0];
 line(dx',dy','linewidth',2) % draw the deformed net
end
fill([-1 0 1 -1],[0 0 -1 -1],'w') % draw the foundation
fill([nos nos+1 nos+1 nos-1],[0 0 -1 -1],'w')
axis off
return
end

function x = gauss(S,f)
n = length(f);
S = [S f]; % Augment S with f
for k=1:n-1 % k counts columns
    [~, rnum] = max(abs(S(k:n,k))); %maximum number given in the quiz
    rnum = rnum + k - 1;
S([rnum k],:) = S([k rnum],:); % Row Swapping
for i = k+1:size(S,2)-1 
    eq1 = S(k,:); %equation 1
    eq2 = S(i,:); %equation 2
    
    m = eq2(k)/eq1(k); %multipl
    
    eq2 = eq2-m*eq1;
    
    S(i,:) = eq2;
end
f = S(:,end);
S = S(:,1:end-1);
end
% Row Swapping
%S(j,:) = S(j, :) + magicnumber*S(k,:);
x = trisolve(S,f);
end
function trisolve(S,f)
x = zeros(n,1);
x(n) = f(n)/S(n,n);
for j = n-1:-1:1
    tmp = 0;
    for k = j+1:n
        tmp = tmp + S(j,k)*x(k);
    end
    x(j) = (f(j)-tmp)/S(j,j);
end
end


