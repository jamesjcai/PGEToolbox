function [a]=trim_gtex_id(a)


for i=1:length(a)
    x=strfind(a{i},'-');
    a{i}=a{i}(1:x(2)-1);
end