function mktestgui(aln)
%MKTESTGUI - Invoke DPRS table for MK test

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[Ds,Ps,Dn,Pn] = dspsdnpn(aln);
	gcbox=dprstable;
	handles = guihandles(gcbox);
	set(handles.ds, 'String', num2str(Ds));
	set(handles.ps, 'String', num2str(Ps));
	set(handles.dr, 'String', num2str(Dn));
	set(handles.pr, 'String', num2str(Pn));
	set(handles.dsps,'String',num2str(Ds+Ps));
	set(handles.drpr,'String',num2str(Dn+Pn));
	set(handles.dsdr,'String',num2str(Ds+Dn));
	set(handles.pspr,'String',num2str(Ps+Pn));
	set(handles.total,'String',num2str(Ds+Ps+Dn+Pn));

