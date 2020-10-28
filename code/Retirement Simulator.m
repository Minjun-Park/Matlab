%Minjun Park, CAAM 210, 11/27/15
%
%Description: With allocated amounts of stocks, bonds, bank, inverted, and
%casino, this project returns savings in the end of each year over 40 years
%
%Example: firstyear = oneyear(7000, 0, 0, 3000, 0)
%
%USAGE: decide how much money to invest in each category(stocks, bonds,
%bank, inverted, casino). considering the probability of type of the year and
%the probability within the category, make strategy
%
function oneyeardriver
for j =1:20 %repeat 20 times
    
n = 1:40; %over 40 years
firstyear = oneyear(7000, 0, 0, 3000, 0); %first year investment
for i = n;  
    if i < 10 %from year 1 to 9
    nextyear = oneyear(firstyear*0.9, firstyear*0.1, 0, 0, 0); %Strategy 1
    firstyear = nextyear+10000; %adding 10000 every year
    end
    if i >=10 && i <15 %from year 10 to 14
    nextyear = oneyear(firstyear*0.7, firstyear*0.1, firstyear*0.2, 0, 0); %Strategy 2
    firstyear = nextyear+10000; %adding 10000 every year
    end
    if i>=15 && i <= 40 %from year 15 to 40
    nextyear = oneyear(firstyear*0.6, firstyear*0.1, firstyear*0.3, 0, 0); %Strategy 3
    firstyear = nextyear+10000; %adding 10000 every year
    end
    %nextyear = 10000 + nextyear
    average(i) = nextyear; % value of end of the year
end

plot(n, average) %years in x-axis, money in y-axis
title('Retirement Savings over 40 years') %title
xlabel('Years') %x-axis
ylabel('Money(in dollars)') %y-axis
hold on
end

c=0; %initialization
a = [];
for j = 1:2000 %repreat 2000 times
n = 1:40;
firstyear = oneyear(7000, 0, 0, 3000, 0); 
for i = n  
    if i < 25
    nextyear = oneyear(firstyear*0.9, firstyear*0.1, 0, 0, 0);
    firstyear = nextyear + 10000;
    end
    if i >=25 && i <35
    nextyear = oneyear(firstyear*0.7, firstyear*0.1, firstyear*0.2, 0, 0);
    firstyear = nextyear + 10000;
    end
    if i>=35 && i <= 40
    nextyear = oneyear(firstyear*0.6, firstyear*0.1, firstyear*0.3, 0, 0);
    firstyear = nextyear + 10000;  
    end
    average(i) = nextyear;
end %{same code as precedent}
k = average(40); %outcome of the 40th year
if k > 1500000 %if the 40th year is greater than 1500000
    c = c+1; %increment c
end
a = [average(40) a]; %getting values of the last years in 2000 repetitions
end
mx = max(a); %maximum
mn = min(a); %minimum
med = median(a); %median
ave = mean(a); %mean
disp('Trials: 2000')
disp('Maximum Ending Value:')
disp(mx)
disp('Median Ending Value:')
disp(med)
disp('Mean Ending Value:')
disp(ave)
disp('Minimum Ending Value:')
disp(mn)
disp('Proportion above 1500000:')
disp(c/2000) % getting proportion that is greater than 1500000
end

function [outputdollars] = oneyear(stocks, bonds, bank, inverted, casino)

yearrand = rand(); % determines how the particular year went
bucketrand = rand(); %getting random number from 0 to 1
if yearrand < 0.05
    yeartype=0; %catastrophic year
elseif yearrand < 0.50
    yeartype =1;%average year
elseif yearrand < 0.95
    yeartype =2; %Good year
elseif yearrand<= 1
    yeartype =3; %Steller year
end

if yeartype == 0 %catastrophic year
    if bucketrand < 0.2 % stocks change
        stocks = stocks*0.6; %reduced by 40%
    elseif bucketrand < 0.6 
        stocks = stocks*0.92; %reduced by 8%
    elseif bucketrand < 0.8
        stocks = stocks*0.96; %reduced by 4%
    elseif bucketrand <= 1
        stocks; %value stays the same
    end
    
    if bucketrand < 0.2 %bonds change
        bonds = bonds*0.6;
    elseif bucketrand < 0.6
        bonds = bonds*0.98;
    elseif bucketrand < 0.8
        bonds = bonds*0.99;
    elseif bucketrand <= 1
        bonds;
    end
    %{same rule applies as before, and the rest has the same rule}
    if bucketrand < 0.2 %inverted change
        inverted = inverted*0.94;
    elseif bucketrand < 0.6
        inverted;
    elseif bucketrand < 0.8
        inverted = inverted*1.06; %increased by 6%
    elseif bucketrand <= 1
        inverted = inverted*1.12; %increased by 12%
    end
end

if yeartype ==1 %Average year
     if bucketrand < 0.2 % stocks change
        stocks = stocks*0.96;
    elseif bucketrand < 0.6
        stocks;
    elseif bucketrand < 0.8
        stocks = stocks*1.04;
    elseif bucketrand <= 1
        stocks = stocks*1.10;
    end
    
    if bucketrand < 0.2 %bonds change
        bonds = bonds*0.98;
    elseif bucketrand < 0.6
        bonds = bonds*1.01;
    elseif bucketrand < 0.8
        bonds = bonds*1.02;
    elseif bucketrand <= 1
        bonds = bonds*1.03;
    end
    
    if bucketrand < 0.2 %inverted change
        inverted = inverted*0.94;
    elseif bucketrand < 0.6
        inverted;
    elseif bucketrand < 0.8
        inverted = inverted*1.06;
    elseif bucketrand <= 1
        inverted = inverted*1.12;
    end    
end

if yeartype ==2 %Good year
     if bucketrand < 0.2 % stocks change
        stocks = stocks*0.96;
    elseif bucketrand < 0.6
        stocks=stocks*1.04;
    elseif bucketrand < 0.8
        stocks = stocks*1.08;
    elseif bucketrand <= 1
        stocks = stocks*1.12;
    end
    
    if bucketrand < 0.2 %bonds change
        bonds = bonds*0.99;
    elseif bucketrand < 0.6
        bonds = bonds*1.02;
    elseif bucketrand < 0.8
        bonds = bonds*1.04;
    elseif bucketrand <= 1
        bonds = bonds*1.06;
    end
    
    if bucketrand < 0.2 %inverted change
        inverted = inverted*0.95;
    elseif bucketrand < 0.6
        inverted=inverted*1.02;
    elseif bucketrand < 0.8
        inverted = inverted*1.03;
    elseif bucketrand <= 1
        inverted = inverted*1.04;
    end    
end

if yeartype ==3 %Steller year
     if bucketrand < 0.2 % stocks change
        stocks = stocks*1.08;
    elseif bucketrand < 0.6
        stocks = stocks*1.16;
    elseif bucketrand < 0.8
        stocks = stocks*1.24;
    elseif bucketrand <= 1
        stocks = stocks*1.50;
    end
    
    if bucketrand < 0.2 %bonds change
        bonds = bonds*1.08;
    elseif bucketrand < 0.6
        bonds = bonds*1.10;
    elseif bucketrand < 0.8
        bonds = bonds*1.12;
    elseif bucketrand <= 1
        bonds = bonds*1.20;
    end
    
    if bucketrand < 0.2 %inverted change
        inverted = inverted*0.96;
    elseif bucketrand < 0.6
        inverted = inverted*1.01;
    elseif bucketrand < 0.8
        inverted = inverted*1.02;
    elseif bucketrand <= 1
        inverted = inverted*1.03;
    end    
end
%{this data is given by the project}
%{each category has own probability that corresponds to the type of the year}

if bucketrand < 0.999 
    bank = bank*1.01;%1 percent increase in saving
else bank = bank*0.8; %loosing 20% with the chance of 0.01%
end
%{bank deposit returns}

if bucketrand < 0.5
    casino = casino*1.90; % getting 90% increase
else casino = 0; %losing all the money
end
%{casino returns}

outputdollars = stocks+bonds+bank+inverted+casino; %this gives the final value of total money after one year, by adding all 4 categories
end
%{this function is used for the strategy to calculate different probabilities}