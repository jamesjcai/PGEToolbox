function [exedir,dlgshown]=pge_getprgmdir(prefstr)

% e.g., prefstr='sweepfinderrun_prgmdir';
exedir='';
dlgshown=false;
dirok=false;
if ispref('pgetoolbox',prefstr)
    exedir=getpref('pgetoolbox',prefstr);
    if isdir(exedir), dirok=true; end
end
if ~dirok
        dlgshown=true;
    	ButtonName=questdlg(sprintf('Do you want to set up external program path?'),...
        			    'External Program First-time Run',...
			            'Yes','No','Cancel','Yes');
        if strcmp(ButtonName,'Yes')
           exedir=uigetdir(pwd);
                    if exedir~=0
                       setpref('pgetoolbox',prefstr,exedir);
                       helpdlg('You are ready to run.','External Program Path Saved')
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