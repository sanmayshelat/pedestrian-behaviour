%% Check combined availability of two or more entities

function [combinedavail] = checkAvailability(agenda)

numMembers = size(agenda,1);
avail = sum(agenda,1);
avail = avail==numMembers;

n = size(avail,2);
f = find([true,diff(avail)~=0,true]);
y = zeros(1,n);
y(f(1:end-1)) = diff(f);
y = cumsum(y(1:n))-(1:n);
combinedavail = y.*avail;

end