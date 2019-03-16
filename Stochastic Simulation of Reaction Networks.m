% Minjun Park, CAAM 210, 10/04/15
%
% This project describes how different equations react with each other
% gilldriver function sets initial parameters such as rtab, nr, tinc, tfin,
% tvec, x, and c, and plots two different plots depending on the usage of
% 7th and 8th reaction equation (calls upon mygill function).
% Usage: input different equations, and then plot the graph with the given
% x-axis and y-axis

function gilldriver
rtab = {[1 -1 2 -1 3 1]%R + Pr -> RPr
    [3 -1 1 1 2 1] %RPr -> R + PR
    [3 -1 4 10 1 1 2 1]%RPr -> 10P + R + Pr
    [4 -2 5 1]%P + P -> D
    [5 -1 4 2]%D -> P + P
    [4 -1]%P ->
    [2 -1 5 -1 6 1]%Pr + D -> I
    [6 -1 2 1 5 1]};%I -> Pr + D

nr=4; %number of times we iternate
tinc = 0.05;
tfin = 100;
tvec = 0:tinc:tfin; %uniform time vector
x=[10 1 0 0 0 0]; %x = [R Pr RPr P D I]

c=[2 1 4 2 0.5 0.05 0 0]; %propensity (without reaction 7 and 8)
X=zeros(nr, length(tvec)); %preallocating
for j=1:nr %lay Gillespie runs along the rows of matrix and average down the columns
    [t,dimer] = mygill(tfin,rtab,x,c); 
    X(j,:)= interp1(t, dimer, tvec); %interpolation
end
avg = mean(X,1); %average
dev = std(X); %standard deviation
plot(tvec,avg) %plot time vs dimer count graph
hold on
plot(tvec,(avg+dev),'r') %plot with std dev added
plot(tvec,(avg-dev),'r') %plot wiht std dev subtracted
xlabel('time')
ylabel('dimer count')
hold off

figure

c=[2 1 4 2 0.5 0.05 0.01 0.05]; %propensity (with reaction 7 and 8)
X2=zeros(nr, length(tvec)); %preallocating
for j=1:nr %lay Gillespie runs along the rows of matrix and average down the columns
    [t,dimer] = mygill(tfin,rtab,x,c);
    X2(j,:)= interp1(t, dimer, tvec); %interpolation
end
avg2 = mean(X2,1); %average
dev2 = std(X2); %standard deviation
plot(tvec,avg2) %plot time vs dimer count graph
hold on
plot(tvec,(avg2+dev2),'r') %plot with std dev added
plot(tvec,(avg2-dev2),'r') %plot wiht std dev subtracted
xlabel('time')
ylabel('dimer count')
hold off
end


% This is mygill function.
% mygill function takes in tfin, rtab, x, and c, and results in [t, dimer].
% Within this function, update function is used.
% Example: mygill(100, rtab, [10 1 0 0 0 0], c=[2 1 4 2 0.5 0.05 0 0])



function [t,dimer] = mygill(tfin, rtab, x, c)
t=[0]; %preallocating t
iter=1; %initial iter
dimer=[0]; %preallocating dimer
while t(iter)<tfin %run as long as t(iter) < tfin
    iter=iter+1; %counter
    a = c .* update(x); %likelihood of each reaction happening
    a0 = sum(a); %add components in a
    r1 = rand; %random number between 0 and 1
    tau = (1/a0)*log(1/r1); %calculate when it is likely to occur
    r2 = rand; %random number between 0 and 1
    cumchance = cumsum(a)/sum(a); %weighted probabilities
    b = find(cumchance > r2,1); %determining which reaction will occur
    x(rtab{b}(1:2:end)) =  x(rtab{b}(1:2:end)) + rtab{b}(2:2:end); %update x according to which reaction occured
    t(iter) = t(iter-1) + tau; %update t(iter) by adding calculated tau to prior t(iter)
    dimer(iter) = x(5); %update the dimer count
end
end

% This is the update function.
% update function updates the remaining and produced reactants after each
% reaction occurs. Each element of h is the number of possible scenarios the
% correponding reaction can occur.
% Example: update([10 1 0 0 0 0])

function h = update(x) %h is an updated probability
h(1) = x(1)*x(2);%reaction 1
h(2) = x(3);%reaction 2
h(3) = x(3);%reaction 3
h(4) = x(4)*(x(4)-1)/2;%reaction 4
h(5) = x(5);%reaction 5
h(6) = x(4);%reaction 6
h(7) = x(2)*x(5);%reaction 7
h(8) = x(6);%reaction 8
end
