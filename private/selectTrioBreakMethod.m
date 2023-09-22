function [methodtype] = selectTrioBreakMethod

prompt = sprintf('SELECT TRIO TYPE:\n1 = HapMap Panel: Yoruba-30-trios\n2 = HapMap Panel: CEPH-30-trios\n3 = Every 3rd Individual\n');
methodtype = input(prompt);


%{
items(1).name = 'TRIO TYPE:';
items(1).default = 1;
items(1).linked = [2 3 4];
items(1).values = [];

items(2).name = '1 = HapMap Panel: Yoruba-30-trios';
items(2).default = 1;
items(2).exclusive = [3 4];
items(2).indent = 1;

items(3).name = '2 = HapMap Panel: CEPH-30-trios';
items(3).default = 0;
items(3).exclusive = [2 4];
items(3).indent = 1;

items(4).name = '3 = Every 3rd Individual';
items(4).default = 0;
items(4).exclusive = [2 3];
items(4).indent = 1;

title = 'Trio Breaking Method';
% msg = sprintf(['Please select sequence type and genetic code']);
out = CSEFlagDialog(items, title);


if ~(isempty(out)),
    if out(1).answer
        if(out(2).answer==1)
        methodtype=1;
        elseif(out(3).answer==1)
        methodtype=2;
        elseif(out(4).answer==1)
        methodtype=3;
    end
else
    methodtype=[];
end
else
methodtype=[];
end
%}