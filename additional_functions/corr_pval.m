% Mitchell Chandler, SIO
% Last updated: 23/04/2021

%Compute the correlation coefficient, r, and p-value between two time
%series, x and y, taking into account the eDOF (Thomson and Emery 2014) of
%each time series.

function [r,pval,eDOF] = corr_pval(x,y)
%remove missing (NaN) values - nb. this assumes missing values are in the same place for each time series 
x = rmmissing(x); 
y = rmmissing(y);

N = length(x); %number of points - nb. function requires x and y be the same size

%find eDOF of each time series
    %eDOFx:
A = xcov(x,'normalized'); %compute autocovariance
AA = A(round(length(A)/2):end); %consider just the positive half of the record
max_integral = max(cumtrapz(AA)); %integral timescale approximated by taking the maximum of the autocorrelation integral
T = 2*max_integral; %integral time scale (multiplied by 2 to take into account +ve and -ve lags)
eDOFx = N/T; 
    %eDOFy:
A = xcov(y,'normalized'); %compute autocovariance
AA = A(round(length(A)/2):end); %consider just the positive half of the record
max_integral = max(cumtrapz(AA)); %integral timescale approximated by taking the maximum of the autocorrelation integral
T = 2*max_integral; %integral time scale (multiplied by 2 to take into account +ve and -ve lags)
eDOFy = N/T;

%effective degrees of freedom taken to be the smaller eDOF, but not larger than the number of data points
eDOF = min([eDOFx eDOFy N]);

%Compute correlation of x and y
r = corr(x,y);

%Compute p-value using the correlation coefficient and eDOF 
%(adapted from Ian Eisenman, and can confirm p-values using
%https://www.danielsoper.com/statcalc/calculator.aspx?id=44) 
tval=abs(r)*sqrt((eDOF-2)/(1-r^2)); %t-value
pval=2*tcdf(-tval,eDOF-2); %two-tailed significance test
end


