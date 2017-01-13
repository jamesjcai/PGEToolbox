function [execmd,exedir,dlgshown]=pge_getexecmd(prefstr)

% e.g., prefstr='sweepfinderrun_dir';
execmd='';
dlgshown=false;
dirok=false;
if ispref('pgetoolbox',prefstr)
    execmd=getpref('pgetoolbox',prefstr);
    if isdir(execmd), dirok=true; end
end
if ~dirok
        dlgshown=true;
    	ButtonName=questdlg('Do you want to set up external program (required for the first-time run)?',...
			    'Setting Up External Program',...
			    'Yes','No','Cancel','Yes');
        if strcmp(ButtonName,'Yes')
           % execmd=uigetfile(pwd);
           if ispc
            [execmd, exedir] = uigetfile('*.exe', 'Pick an External Program');
           else
            [execmd, exedir] = uigetfile('*.*', 'Pick an External Program');
           end
                    if isequal(execmd,0) || isequal(exedir,0)
                       execmdfull=fullfile(pathname, filename);
                       setpref('pgetoolbox',prefstr,execmdfull);
                       helpdlg('You are ready to run the external program.','External Program Setup Ready')
                    end
        end
        return;
end


%{                                
[selectedButton,dlgShown]=uigetpref('pgetoolbox',... % Group
       'sweepfinderrun_ask',...                      % Preference
       'Saving Directory',...              % Window title
       {'Do you want to save the choosen directory?'},...
       {'always','never';'Yes','No'},...       % Values and button strings
       'ExtraOptions','Cancel',...             % Additional button
       'DefaultButton','Cancel');
switch selectedButton
    case {'always','Yes'}
        setpref('pgetoolbox','sweepfinderrun_dir',workingdir);
        helpdlg('You are ready to run the program.','Executing directory saved')
    case {'never','No','Cancel'}
        % do nothing
end
%}