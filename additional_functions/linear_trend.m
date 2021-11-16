% Mitchell Chandler, SIO
% Last updated: 02/03/2021

%Compute linear trend and assess the confidence interval of the gradient of
%the trend line taking into account any autocorrelation following Thomson
%and Emery (2014).

function [p,Yhat,CI] = linear_trend(X,Y,dt,alpha)
%make sure both X and Y are [N x 1]
if size(Y,2) > 1
    Y=Y';
end
if size(X,2) > 1
    X=X';
end

N = length(Y); %number of points

%linear fit
p = polyfit(X,Y,1); %p(1) is gradient p(2) is intercept
Yhat = polyval(p,X); %linear fit

%compute autocovariance
A = xcov(Y,'normalized');
AA = A(round(length(A)/2):end); %consider just the positive half of the record

%integral timescale approximated by taking the maximum of the autocorrelation integral
[max_integral,idx_integral] = max(cumtrapz(AA));
T = 2*max_integral*dt; %integral time scale (multiplied by 2 to take into account +ve and -ve lags)

%effective degrees of freedom
eDOF = (N*dt)/T;

%maximum DOF is the number of data points
if eDOF > N 
    eDOF = N;
end

%compute CI
s_e = sqrt(1/(N-2)*sum((Y-Yhat).^2));
N_s_x = sqrt(sum((X-mean(X)).^2));
tval = tinv(1-alpha/2,eDOF-2);
CI = (s_e * tval) / (N_s_x);
end

