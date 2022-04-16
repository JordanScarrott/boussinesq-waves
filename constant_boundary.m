
% Sets k many rows and columns at the borders of n to the same as k+1
% within the border

% All +'s will be the same as x when k==2
% ++++++++++++++++++++
% ++++++++++++++++++++
% ++xxxxxxxxxxxxxxxx++
% ++x--------------x++
% ++x--------------x++
% ++x--------------x++
% ++x--------------x++
% ++xxxxxxxxxxxxxxxx++
% ++++++++++++++++++++
% ++++++++++++++++++++

function n = constant_boundary(n, k)
    for i=k:-1:1
        n(:,i) = n(:,i+1);
        n(:,end-(i-1)) = n(:,end-i);
        n(i,:) = n(i+1,:);
        n(end-(i-1),:) = n(end-i,:);
    end
end

