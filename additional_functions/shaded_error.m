% Mitchell Chandler
% Last updated: 08/09/2020

%Plot a line and a shaded error bar around the line.

function [] = shaded_error(x,y,err,colour,a,lwidth)
%default alpha (opacity) and linewidth values if none are given
if ~exist('a','var') %
      a = 0.2;
 end
if ~exist('lwidth','var')
      lwidth = 2;
end
%check if row or column vector and convert to column vector
if isrow(x)
    x=x';
end
if isrow(y)
    y=y';
end
if isrow(err)
    err=err';
end
%plot
hold on
fill([x;flipud(x)],[y-err;flipud(y+err)],colour,'linestyle','none','facealpha',a,'HandleVisibility','off'); %shaded error
plot(x,y,'Color',colour,'LineWidth',lwidth) %line
end

