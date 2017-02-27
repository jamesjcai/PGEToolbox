function [pw1,pw0]=cdpge(isconfirmed)
%CDPGE - Changes to PGEToolbox directory

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

if nargin < 1
    isconfirmed=true;
end
pw0=pwd;
pw1=fileparts(which(mfilename));
if ~strcmp(pw0,pw1) && ~isconfirmed
    [selectedButton]=uigetpref('PGEToolbox',... % Group
           'cdpge_ask',...                               % Preference
           'Changing Working Directory',...              % Window title
           {'Do you want to change current working directory to PGEToolbox directory?'},...
           {'always','never';'Yes','No'},...       % Values and button strings
           'ExtraOptions','Cancel',...             % Additional button
           'DefaultButton','Yes');
    switch selectedButton
        case {'always','Yes'}
            cd(pw1);
        case {'never','No','Cancel'}
            % do nothing
    end
else
    cd(pw1);
end