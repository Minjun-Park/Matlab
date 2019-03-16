%Minjun Park, CAAM 210, 11/15/15
%
%Description: Forms gene networks, which shows how each gene interact with one another 
%
%Usage: with given rules, we will express the action at each gene, aka
%logic table
%
function genestmdriver
s = [1 0 1 0 0 1]; %current state
wire = [4 2 3; 5 3 2; 3 6 1; 5 4 6; 6 1 2; 3 5 6]; %represents each node
n = size(wire,1); %wire-rows, 1-column
a = zeros(n+1); %preallocation of a
for i=1:n 
    a(wire(i,:),i) = 1; 
    ids{i} = num2str(i);
end
ids{n+1} = 'Gene Net'; %title
g = biograph(a,ids); %biograph command
selfcon = find(diag(a)==1); 
for i=1:length(selfcon) 
    g.nodes(selfcon(i)).Shape = 'diamond'; %changing to diamond form 
end
view(g)
%{this is the given code from notebook, which shows how 6 genes interact with each other}

wire = [4 2 3; 5 3 2; 3 6 1; 5 4 6; 6 1 2; 3 5 6];
n = size(wire,1); %preallocating n 
rule = [231; 90; 5; 108; 61; 62]; %node 2 is changed
STM = genestm(wire,rule); %call STM
for i=1:2^n 
    ids{i} = num2str(i-1); %this gives the number of ids we are going to create
end
ids{2^n+1} = 'State Transition Diagram'; %title of the diagram
g = biograph(STM,ids); %biograph gives visualized network
selfcon = find(diag(STM)==1); 
for i=1:length(selfcon) 
    g.nodes(selfcon(i)).Shape = 'diamond'; %arrow pointing to itself makes gene into a diamond shape
end
view(g)
end
%{this is a driver function that creates two distinct images}

function STM = genestm(wire,rule)
STM = zeros(65,65); %preallocating STM
ruletable = zeros(6,8); %we have 6 genes with 2^3 columns
for i = 1:6
   ruletable(i,:) = d2b(rule(i), 8); %creating ruletable in binary form
end
n = size(wire, 1); %preallocating n
for i = 1:2^n
    s= d2b(i-1,n); %s is a binary form
    for k = 1:n
        m = s(wire(k,:)); %calling whole row
        q = b2d(m);%change to decimal
        ns(k) = ruletable(k,8-q); %next stage
    end
    j = b2d(ns); %get decimals of next stage
    STM(i,j+1) = 1; 
end
end 
function b = d2b(r,c) 
b = zeros(1,c);
while r>0
    b(c-floor(log2(r))) = 1;
    r= r-2^(floor(log2(r))); %get each decimals
end
end
%{changing from decimal to binary}
function d = b2d(b)
d = 0;
for i = 1:length(b)
    d = d + b(i)*2^(length(b) - i); %covert to either 0 or 1
end
end
%{changing from binary to decimal}


