function [y] = chkmbeinstall(required)
if nargin < 1
    required = true;
end
y = exist('cdmbe.m', 'file') ~= 0;
if required && ~y

    disp('This function requires MBEToolbox.');
        disp('Using the code below to install MBEToolbox:');
        disp(' ');
        disp('unzip(''https://github.com/jamesjcai/MBEToolbox/archive/master'');');
        disp('addpath(''./MBEToolbox-master'');');
        disp(' ');
            error('This function requires MBEToolbox.');
        end