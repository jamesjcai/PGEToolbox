function [tag] = sigtag(p)
%SIGTAG - 

if (p>0.01 && p<0.05)
    tag=' (*significant, 0.01 < P < 0.05)';
elseif (p>0.001 && p<=0.01)
    tag=' (**significant, 0.001 < P < 0.01)';
elseif (p<=0.001)
    tag=' (***significant, P < 0.001)';
elseif (p>0.05 && p<=0.1)
    tag=' (not significant, P > 0.05)';
else
    tag=' (not significant, P > 0.1)';
end
