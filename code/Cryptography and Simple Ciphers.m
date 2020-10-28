% Minjun Park, CAAM 210, 10/14/15
%
% This is decoder function
% This code calls encrypted file, and using letter transition
% probabilities, it decodes the file. The logic of the function is simply
% replacing a pair of letters with the higher probability 
% 
% USAGE: Call three different filename, and then this function converts
% those mysterious document into sentences that make sense
function decoderdrive 
code1 = fileread('encodedtext1.txt'); % calls the message
code2 = fileread('encodedtext2.txt'); % message 2
code3 = fileread('encodedtext3.txt'); % message 3
decoder(code1) %run decoder function 
decoder(code2) %run seoncd one
decoder(code3) % run the third one
end


function e = decoder(x) 
a = downlow(x); %calles the function
y = randperm(27); %randomly display 27 numbers
for i=1:10^5 %enough iteration to check the probability
    ymaybe = y; %replaces the ymaybe with y
    k = y(r); %
    r = randperm(numel(y),2); %randomly select one pair
    l = k(1); %if it is equal to the first element
    k(2) = l; %switching the pair
    k(1) = k(2); %repalce first element with the second one
    ymaybe(r(1))=k(1); 
    ymaybe(r(2))=k(2); 
     if loglike(y,a) < loglike(ymaybe,a) %when the probability of ymaybe is bigger
          y = ymaybe; %replace the y with ymaybe
           elseif loglike(y,a) > loglike(ymaybe,a)  %when the probability of y is bigger
               if rand < exp(-loglike(y,a)+loglike(ymaybe,a)) %this is the given formula by the assignment
                    y = ymaybe; %ymaybe becomes y
               end
     end
end
for i=1:numel(a)
       a(i)=y(a(i));      
end
e = downlowinv(a);
end
%call decoder three times, and then print the output to the screen.

function o = downlow(x)
o = zeros(1, numel(x));
for i = 1:numel(x)
    if x(i) == 96
        o(i) = double(x(i)) - 69; %subtract to account for space
    else o(i) = double(x(i)) - 96; %subtract to get regular alphabets
    end
end
end
% this function converts each of the ascii codes corresponding to the
% letters a through z and space to the numbers 1 through 27



function v = downlowinv(x)
v = zeros(1, numel(x)); %set array
for i = 1:numel(x);
   if x(i) ==  27 %this is the space
       v(i) = 32; %then this gives ' ' for switching
   else
       v(i) = x(i) + 96; %this gives regular alphabet
   end
end
v = char(v);
end
% this functions takes a number 1 to 27 and returns the corresponding
% character



function k = loglike(x, z)
k=0; %initialization
j = textread('letterprob.mat'); %calls the this probability function
for i = 2:numel(z) 
    k = k + log(j(x(z(i-1), x(z(i))))); %this calls the probablity of each letter showing up after another letter
end
end
% this function takes in the coded text(as an array of numbers 1 through
% 27) and my guess(a 1*27 array) which computes the log-likelihood of my
% guess given the inputted text
%
%
%
%
%
%



