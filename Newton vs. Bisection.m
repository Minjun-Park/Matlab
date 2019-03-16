%Minjun Park, CAAM 210, 09/21/15
%
%This function shows comparionson between Bisection and Newton's method
%
%USUAGE: X-axis would be tolerance and Y-axis is # of iterations, and it
%shows what happens as tolerance goes smaller
function solverace
L = 1;
a = 0.1;
b = 3;
x = (a + b)/2;
 
for j = 1:8
    tol(j) = 10^(-j); %This gives tolerance of the x-axis
    [var, iter(j)] = denewt(x,tol(j),L); %takes the variables
    [~, iteration(j)] = debis(a,b,tol(j),L); %takes the variables
end
semilogx (tol, iter, 'o-') %this is for denewt function
hold on
semilogx (tol, iteration, 'x-') %this is for debis function
title('Bisection vs. Newton') 
xlabel('Tolearnce') %shows name of x-axis
ylabel('# of iterations') %shows name of y-axis 
legend('Newton', 'Bisection') %shows legend the plot
grid on
end

function [x, iter] = denewt(x,tol,L)
    %This function fcodes Newton's method
    iter = 0;
   while abs(coolfun(x,L)) > tol 
       x = x - coolfun(x, L)/coolfundx(x, L);
       iter = iter + 1;
   end
end

function [x,iter] = debis(a,b,tol,L)

x = NaN;                % not a number
iter = 0;               % Initialize iteration counter

if coolfun(a,L)*coolfun(b,L) > 0,   % check to see if [a,b] contains a root
   disp('Sorry, I can not be sure there is a root');
end
x = (a+b)/2;                    % make first midpoint
while abs(coolfun(x,L)) > tol 
      if coolfun(a,L)*coolfun(x,L) < 0      % if f crosses zero between a and x then
         b = x;                 % move b
      else
         a = x;                 % move a
      end
      x = (a+b)/2;              % make new midpoint
      iter = iter + 1;          % Update iteration counter
end   % while
end

function val = coolfun(x,L)
val = sin(x*L) + x*cos(x*L);
end
%{this function returns cooling bar function}

function val = coolfundx(x,L) 
    %Evaluate the derivetive, with respect to x, of coolfun
val = L*cos(x*L) + cos(x*L) - x*L*sin(x*L);
end
%{this function takes derivative of coolfun, and shows the equation}
