% Mitchell Chandler, SIO
% Last updated: 03/03/2021

%Filter the data X using the window W [e.g. W =
%(1/sum(triang(w_size)))*triang(w_size)] with a window of length w_size.

function [X_filtered] = conv_filt(X,W,w_size)
%apply filter
filtered_data = conv(W,X); 
%cut to size of original data
X_filtered = filtered_data(ceil(w_size/2):end-floor(w_size/2)); 
%set values where the full window cannot be applied as NaNs
X_filtered(1:floor(w_size/2)) = NaN; 
X_filtered((end+1)-floor(w_size/2):end) = NaN;
end