function [b] = mat2cellstr(a)
%http://www.mathworks.fr/matlabcentral/newsreader/view_thread/156758

%b = strread(num2str(a),'%s');
b=textscan(num2str(a),'%s');
b=b{1};