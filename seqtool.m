function varargout = seqtool(varargin)
%seqtool - PGEToolbox GUI
%      seqtool, by itself, creates a new seqtool or raises the existing
%      singleton*.
%
%      H = seqtool returns the handle to a new seqtool or the handle to
%      the existing singleton*.
%
%      seqtool('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in seqtool.M with the given input arguments.
%
%      seqtool('Property','Value',...) creates a new seqtool or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PGEGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seqtool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $
%

% Edit the above text to modify the response to help seqtool

% Last Modified by GUIDE v2.5 29-Dec-2012 18:22:29
global pgeversionstr;
pgeversionstr='3.1';


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seqtool_OpeningFcn, ...
                   'gui_OutputFcn',  @seqtool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before seqtool is made visible.
function seqtool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seqtool (see VARARGIN)

% Choose default command line output for seqtool
global pgeguiolddir;
handles.output = hObject;

scrsz=get(0,'ScreenSize');
pos_act=get(gcf,'Position');
xr=scrsz(3)-pos_act(3);
xp=round(xr/2);
yr=scrsz(4)-pos_act(4);
yp=round(yr/2);
set(gcf,'position',[xp yp pos_act(3) pos_act(4)]);

if ~ispref('PGEToolbox','lastworkingdir')
   addpref('PGEToolbox','lastworkingdir',pwd)
else
   try
    cd(getpref('PGEToolbox','lastworkingdir',pwd));
   catch ME
    rmpref('PGEToolbox','lastworkingdir');
   end

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes seqtool wait for user response (see UIRESUME)
% uiwait(handles.PGEToolboxGUI);
% load 'seqtool.mat' PGEGUI_SaveDistance
%dirstr = fileparts(which(mfilename));
%cd(dirstr);
pgeguiolddir=pwd;
cdpge;
%warning off all
% deletetempfiles;
%clear global aln aln_ori
SetMenuStatus(handles);


% --- Outputs from this function are returned to the command line.
function varargout = seqtool_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global aln;
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Data_Callback(hObject, eventdata, handles)
% hObject    handle to Sequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Sequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFASTAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
aln2=pge_openfile;
if ~(isempty(aln2)),
    aln=aln2;
    if (hasgap(aln))
        aswer=questdlg('Do you want to remove gaps (completely delete)?');
     switch aswer,
     case 'Yes',
           aln_ori=aln;
           aln=rmgaps(aln);
     case 'No',

     case 'Cancel',
          aln=[];
      otherwise
          aln=[];
     end
    end
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function SaveFASTAFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFASTAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
writefasta(aln);

% --------------------------------------------------------------------
function SavePhylipFile_Callback(hObject, eventdata, handles)
% hObject    handle to SavePhylipFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
writephylip_i(aln);


% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pgeguiolddir aln
aln=[];
cd(pgeguiolddir);
clear global aln olddir pgeversionstr
warning on all
close;

% --------------------------------------------------------------------
function CloseDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAln (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aln=[];
clear global aln
SetMenuStatus(handles);

% --------------------------------------------------------------------
function ViewSequences_Callback(hObject, eventdata, handles)
% hObject    handle to ViewSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
viewseq(aln);

% --------------------------------------------------------------------
function ViewSequencesInfo_Callback(hObject, eventdata, handles)
% hObject    handle to ViewSequencesInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
aln
try
   pge_viewdata(aln)
catch ME
   errordlg(ME)
end

% --------------------------------------------------------------------
function RemoveGapsInSequences_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveGapsInSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
if (hasGap(aln))
    aln_ori=aln;
    aln=rmgaps(aln); viewseq(aln);
    SetMenuStatus(handles);
else
    return;
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function IncludeExcludeSequences_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
[s,v] = choosebox('Name','Pick outgroup sequence(s)','PromptString','Sequences available:',...
    'SelectString','Selected sequences:','ListString',aln.seqnames);
if (v==1)
    aln_ori=aln;
    aln.seq = aln.seq(s,:);
    aln.seqnames=aln.seqnames(s);
    SetMenuStatus(handles);
else
    return;
end
end

%{
 --------------------------------------------------------------------
function WorkOnReverseStrand_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnReverseStrand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=revseq(aln);
disp('Working on revserse of sequences now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end
%}

%{
 --------------------------------------------------------------------
function WorkOnRevcomStrand_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnRevcomStrand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=revcomseq(aln);
disp('Working on reversed complement of sequences now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end
%}



% --------------------------------------------------------------------
function WorkOnPloymorphicSites_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnPloymorphicSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=extractsegregatingsites(aln);
if isempty(aln.seq)
    warndlg('No segregating site.')
    aln=aln_ori;
    aln_ori=[];
else
    disp('Working on ploymorphic (segregating) sites only!');
end
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function WorkOnParsimonyInformativeSites_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnParsimonyInformativeSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
[aln,informativesites]=extractinformativesites(aln);
if isempty(informativesites)
    warndlg('No informative site.')
    aln=aln_ori;
    aln_ori=[];
else
    disp('Working on Parsimony informative sites only!');
end
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function RestoreSequences_Callback(hObject, eventdata, handles)
% hObject    handle to RestoreSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln=aln_ori;
disp('Sequences restored!');
aln_ori=[];
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ReportPolymorphicSites_Callback(hObject, eventdata, handles)
% hObject    handle to ReportPolymorphicSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
   reportpolysites(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Tajima89Test_Callback(hObject, eventdata, handles)
% hObject    handle to Tajima89Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
   tajima89d_test(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function HelpContents_Callback(hObject, eventdata, handles)
% hObject    handle to HelpContents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpwin pgetoolbox;


% --------------------------------------------------------------------
function Demos_Callback(hObject, eventdata, handles)
% hObject    handle to Demo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
demo pgetoolbox

% --------------------------------------------------------------------
function SaveMatFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
    [filename, pathname, filterindex] = uiputfile( ...
       {'*.mat', 'MAT-file (*.mat)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];

	if (filterindex==1)
	if (isempty(find(filename=='.', 1)))
		filename=[filename,'.mat'];
	end
	end
save(filename,'aln');


%{
 --------------------------------------------------------------------
function LoadMatFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
    % try: uiopen('load')    

    [filename, pathname] = uigetfile( ...
       {'*.mat', 'MAT-file (*.mat)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Mat file');
	if ~(filename), return; end
	filename=[pathname,filename];
alninmat = load(filename);

if ~(isfield(alninmat,'aln'))
    warndlg ('Invalid MAT file for alignment input')
    return;
else
    aln=alninmat.aln;
    SetMenuStatus(handles);
end
%}

% --------------------------------------------------------------------
function HelpWebSite_Callback(hObject, eventdata, handles)
% hObject    handle to HelpWebSite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    web('http://bioinformatics.org/pgetoolbox/','-browser')
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function McDKTestCMD_Callback(hObject, eventdata, handles)
% hObject    handle to McDKTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    mktestcmd(aln);
catch ME
    errordlg(ME.message)
end

%if (length(unique(aln.population))==2)
%    taxset1=find(aln.population==1);
%    taxset2=find(aln.population==2);
%else
%    [s,v] = choosebox('Name','Interspecies Sequences','PromptString',...
%    'Sequences from species 1:','SelectString','Sequences from species 2:',...
%    'ListString',aln.seqnames);
%    if (v==1)
%    taxset2=s';
%    [n,m]=size(aln.seq);
%    taxset1=setdiff([1:n],taxset2);
%    mktestcmd(aln,taxset1,taxset2);
%    % SetMenuStatus(handles);
%    else
%        return;
%    end;
%end
%mktestcmd(aln,taxset1,taxset2);


% --------------------------------------------------------------------
function CheckLatestVersion_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLatestVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pgeversionstr

try
	latestver = urlread('http://www.bioinformatics.org/pgetoolbox/LatestVersionTag');
catch
	warndlg('Can''t read from URL.');
	return;
end

if (strcmp(latestver,pgeversionstr)==1),
    helpdlg('You are running the latest version.','No update available')
else
	ButtonName=questdlg(sprintf(['Download latest version %s?'], latestver), ...
			    'Check Latest Version', ...
			    'Go to website','Cancel','Go to website');
	switch ButtonName,
	    case 'Go to website',
           web('http://bioinformatics.org/pgetoolbox/','-browser')
	    otherwise
		return;
	end
end



%{
 --------------------------------------------------------------------
 function CiteSoftware_Callback(hObject, eventdata, handles)
% hObject    handle to CiteSoftware (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	ButtonName=questdlg('Looking at the reference?', ...
			    'Cite PGEToolbox in your research', ...
			    'Go to website','Cancel','Go to website');
	switch ButtonName,
	    case 'Go to website',
           web('http://www.biomedcentral.com/1471-2105/6/64','-browser')
	    otherwise
		return;
	end
%}

% --------------------------------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function GenerateRandomSequences_Callback(hObject, eventdata, handles)
% hObject    handle to GenerateRandomSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
   prompt={'Number of Sequences:','Length of Sequences:'};
   def={'10','300'};
   dlgTitle='Generate Random Sequences';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
      n=str2num(answer{1});
      m=str2num(answer{2});
      try
          sroot=randseq(1,m);
          s=repmat(sroot,n,1);
	  seedms_mex;
          [segseq,segsite]=ms_mex(n,0,round(m/2));
          segsite=ceil(segsite*m);
          alle2=randseq(1,m);
          for k=1:size(segseq,2)
            colk=uint8(segseq(:,k));
            colk(colk==1)=alle2(k);
            colk(colk==0)=sroot(k);
            s(:,segsite(k))=colk;
          end
      catch exception
          disp(exception.identifier);
          s=randseq(n,m);
      end
      names={};

   for k=1:n,
       names{k}=['seq',sprintf('%d',k)];
   end


   aln.seqtype = 1;
   %aln.geneticcode = geneticcode;
   aln.seqnames = names;
   aln.seq=s;
   	aln.locus=ones(n,1);
	aln.population=ones(n,1);
	aln.count=ones(n,1);
   SetMenuStatus(handles);
   aln
end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function AssignSequencesInWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to AssignSequencesInWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    %{
    if ~isempty(aln),
    assignin('base','aln',aln);
    %disp('Sequences have been assigned as the variable, aln, in workspace.')
    helpdlg('ALN assigned.','Assign varibles into workspace')
    else
    warndlg('Sequences is empty.','Assign varibles into workspace')
    end
    %}
    export2wsdlg({'Save alignment (ALN) to variable named:'},...
                 {'aln'},{aln},...
                 'Export to Workspace');
catch ME
  errordlg(ME)
end


% --------------------------------------------------------------------
function SequenceTextEditor_Callback(hObject, eventdata, handles)
% hObject    handle to SequenceTextEditor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
   prompt={'FASTA format:'};
   def={''};
   dlgTitle='Sequences Text Editor';
   lineNo=6;
   AddOpts.Resize='off';
   AddOpts.WindowStyle='modal';
   AddOpts.Interpreter='none';
   answer=inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);

 if ~(isempty(answer)),
      seq=cellstr(answer{1});
      n=length(seq);
 if (n>=2&&~isempty(seq(2)))
       tempf = tempname;
       [fid,Msg] = fopen(tempf,'wt');
       if fid == -1, error(Msg); end
       for k=1:n
           fprintf(fid,'%s\n',char(seq(k)));
       end
       fclose(fid);
   aln=readfasta(tempf);
   SetMenuStatus(handles);

    n=size(aln.seq,1);

    aln.locus=ones(n,1);
	aln.population=ones(n,1);
	aln.count=ones(n,1);
    aln
else
    helpdlg('Sorry, not valid sequences.')
end
end

catch ME
    errordlg(ME.message)
end



function SetMenuStatus(handles)
global aln aln_ori

%set(handles.McDKTest,'enable','off');
set(handles.HKATest,'enable','off');
set(handles.MKRPFTest,'enable','off');

proteinonly=[];
nucleotideonly=[];

alltags=[handles.AboutPGEToolbox,...
handles.Analysis,...
handles.AssignLocus,...
handles.AssignPopulation,...
handles.AssignSequencesInWorkspace,...
handles.CheckLatestVersion,...
handles.Citation,...
handles.CloseDataFile,...
handles.CoalescentSimulations,...
handles.Data,...
handles.Demos,...
handles.EstimateTheta,...
handles.Exit,...
handles.FayWusTestH,...
handles.FayWycoffWu01,...
handles.File,...
handles.FoldedSiteFrequencySpectrum,...
handles.FourGameteTest,...
handles.FuFsTest,...
handles.FuLi93Tests,...
handles.GenerateRandomSequences,...
handles.HaplotypeDiversity,...
handles.Help,...
handles.HelpContents,...
handles.HelpWebSite,...
handles.IncludeExcludeSequences,...
handles.IncludeExcludeSites,...
handles.LDDistancePlot,...
handles.LDDistancePlotD,...
handles.LDDistancePlotDAbs,...
handles.LDDistancePlotDPrime,...
handles.LDDistancePlotDPrimeAbs,...
handles.LDDistancePlotR,...
handles.LDDistancePlotR2,...
handles.LDMatrixPlot,...
handles.LDMatrixPlotD,...
handles.LDMatrixPlotDAbs,...
handles.LDMatrixPlotDPrime,...
handles.LDMatrixPlotDPrimeAbs,...
handles.LDMatrixPlotR,...
handles.LDMatrixPlotR2,...
handles.LDOverview,...
handles.LinkageDisequilibrium,...
handles.McDKTestCMD,...
handles.McDKTestGUI,...
handles.MKRPFTest,...
handles.MKRPFTestExampleInput,...
handles.MKRPFTestWriteInputLine,...
handles.NucleotideDiversity,...
handles.OpenDataFile,...
handles.R2Test,...
handles.RatioOfAdaptiveNonsynSub,...
handles.RatioOfAdaptiveNonsynSubFWW,...
handles.RatioOfAdaptiveNonsynSubSEW,...
handles.Recombination,...
handles.RemoveGapsInSequences,...
handles.ReportPolymorphicSites,...
handles.RestoreSequences,...
handles.SaveExportDataAs,...
handles.SaveFASTAFile,...
handles.SaveMatFile,...
handles.SavePhylipFile,...
handles.SequenceTextEditor,...
handles.SiteFreqSpectrum,...
handles.SlidingWindowsAnalysis,...
handles.SlidingWinFuDStar,...
handles.SlidingWinFuFs,...
handles.SlidingWinTajimaD,...
handles.SmithEyreWalker02,...
handles.SNPTool,...
handles.Tajima89Test,...
handles.TestOfIndependence2x2Table,...
handles.Tools,...
handles.UnfoldedSiteFrequencySpectrum,...
handles.ViewSequences,...
handles.ViewSequencesInfo,...
handles.WallBQ,...
handles.Watterson78F,...
handles.WorkOnParsimonyInformativeSites,...
handles.WorkOnPloymorphicSites];



cannotemptyaln = [
    handles.SaveExportDataAs,...
    handles.CloseDataFile,...
    handles.MKRPFTestWriteInputLine,...
    handles.FayWycoffWu01,...
    handles.SmithEyreWalker02,...
    handles.ViewSequences,...
    handles.ViewSequencesInfo,...
    handles.AssignLocus,...
    handles.AssignPopulation,...
    handles.RestoreSequences,...
    handles.RemoveGapsInSequences,...
    handles.WorkOnParsimonyInformativeSites,...
    handles.WorkOnPloymorphicSites,...
    handles.IncludeExcludeSequences,...
    handles.IncludeExcludeSites,...
    handles.ReportPolymorphicSites,...
    handles.NucleotideDiversity,...
    handles.EstimateTheta,...
    handles.SiteFreqSpectrum,...
    handles.Tajima89Test,...
    handles.FayWusTestH,...
    handles.FuLi93Tests,...
    handles.SlidingWindowsAnalysis,...
    handles.HaplotypeDiversity,...
    handles.Watterson78F,...
    handles.FuFsTest,...
    handles.R2Test,...
    handles.McDKTestCMD,...
    handles.McDKTestGUI,...
    handles.Recombination,...
    handles.LinkageDisequilibrium,...
    handles.WallBQ];

if (isempty(aln)),
    set(cannotemptyaln,'enable','off');
else
   if (hasgap(aln))
        set(cannotemptyaln,'enable','on');
        set(handles.Analysis,'enable','off');
        set(handles.FayWycoffWu01,'enable','off');
        set(handles.SmithEyreWalker02,'enable','off');
   else
        set(cannotemptyaln,'enable','on');
   end

   if ~(isempty(aln_ori)),
        set(handles.RestoreSequences,'enable','on');
    else
        set(handles.RestoreSequences,'enable','off');
   end
end
% Update handles structure
% guidata(hObject, handles);


% --------------------------------------------------------------------
function FuLi93Tests_Callback(hObject, eventdata, handles)
% hObject    handle to FuLi93DS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
fuli93dsfs_test(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function EstimateTheta_Callback(hObject, eventdata, handles)
% hObject    handle to EstimateTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
   estimatetheta(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function McDKTestGUI_Callback(hObject, eventdata, handles)
% hObject    handle to DPRSTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    mktestgui(aln)
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function NucleotideDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to NucleotideDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    nucdiv(aln);
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function CoalescentSimulations_Callback(hObject, eventdata, handles)
% hObject    handle to CoalescentSimulations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
coalsimdlg
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Citation_Callback(hObject, eventdata, handles)
% hObject    handle to Citation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(' ')
disp('Cai JJ')
disp('PGEToolbox: a MATLAB toolbox for population genetics and evolution.')
disp('http://bioinformatics.org/pgetoolbox/')
disp(' ')

% --------------------------------------------------------------------
function AboutPGEToolbox_Callback(hObject, eventdata, handles)
% hObject    handle to AboutPGEToolbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pgeversionstr
info{1}='';
info{2}='PGEToolbox - Population Genetics & Evolution Toolbox';
info{3}='';
info{4}=['Version ',pgeversionstr];
info{5}='';
info{6}='Author: James Cai';
info{7}='';
info{8}='Copyright 2010-2013  All Rights Reserved';
info{9}='';
info{10}='Please reference this software when using as part of';
info{11}='research.';
info{12} ='';
helpdlg(info,'About');


% --------------------------------------------------------------------
function Watterson78F_Callback(hObject, eventdata, handles)
% hObject    handle to Watterson78W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    watterson78f(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function R2Test_Callback(hObject, eventdata, handles)
% hObject    handle to R2Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    r2_test(aln);
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function WallBQ_Callback(hObject, eventdata, handles)
% hObject    handle to WallBQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
try
    wall99bq(aln);
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function FuFsTest_Callback(hObject, eventdata, handles)
% hObject    handle to FuFsTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
try
    fu97fs(aln);
catch ME
    errordlg(ME.message)
end



%{
 --------------------------------------------------------------------
function StrobeckSTest_Callback(hObject, eventdata, handles)
% hObject    handle to StrobeckSTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
disp('Using the same function as Fu''s Fs test.')
try
    Fu97Fs(aln);
catch ME
    errordlg(ME.message)
end
%}



% --------------------------------------------------------------------
function SaveExportDataAs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveExportDataAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function SlidingWindowsAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to SlidingWindowsAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [winsiz,stpsiz] = i_getwindowsizestep()
   prompt={'Window size:','Step size:'};
   def={'30','10'};
   dlgTitle='Window Size and Step Size';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 winsiz=0; stpsiz=0;
 if ~(isempty(answer)),
   winsiz=str2num(answer{1});
   stpsiz=str2num(answer{2});
 end


% --------------------------------------------------------------------
function SlidingWinTajimaD_Callback(hObject, eventdata, handles)
% hObject    handle to SlidingWinTajimaD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln

[winsiz,stpsiz] = i_getwindowsizestep;
if ~(winsiz>0 && stpsiz>0)
   return;
else
	try
	    h=figure;
	    slidingfun(@tajima89d_test,aln,winsiz,stpsiz);
        ylabel('Tajima''s D')
    catch
        close(h);
        errordlg(lasterr)
	    %warndlg('Inappropriate WINSIZE or STEP.')
	end
end



% --------------------------------------------------------------------
function SlidingWinFuFs_Callback(hObject, eventdata, handles)
% hObject    handle to SlidingWinFuFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[winsiz,stpsiz] = i_getwindowsizestep;
if ~(winsiz>0 && stpsiz>0)
   return;
else
	try
	    h=figure;
	    slidingfun(@fu97fs,aln,winsiz,stpsiz);
        ylabel('Fu''s Fs')
    catch
        close(h);
        errordlg(lasterr)
	    %warndlg('Inappropriate WINSIZE or STEP.')
	end
end


% --------------------------------------------------------------------
function SlidingWinFuDStar_Callback(hObject, eventdata, handles)
% hObject    handle to SlidingWinFuDStar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln

[winsiz,stpsiz] = i_getwindowsizestep;
if ~(winsiz>0 & stpsiz>0)
   return;
else
	try
	    h=figure;
	    subplot(2,1,1)
	    switchoutput=0;
	    slidingfun(@fuli93dsfs_test,aln,winsiz,stpsiz,switchoutput);
	    ylabel('Fu and Li''s D*')
	    subplot(2,1,2)
	    switchoutput=1;
	    slidingfun(@fuli93dsfs_test,aln,winsiz,stpsiz,switchoutput);
	    ylabel('Fu and Li''s F*')
    catch
        close(h);
        errordlg(lasterr)
	    %warndlg('Inappropriate WINSIZE or STEP.')
	end
end

% --------------------------------------------------------------------
function LinkageDisequilibrium_Callback(hObject, eventdata, handles)
% hObject    handle to LinkageDisequilibrium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LDOverview_Callback(hObject, eventdata, handles)
% hObject    handle to LDOverview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    warning off MATLAB:polyfit:PolyNotUnique
    linkdisequ(aln,0,1);
catch
    warndlg('No pairwise comparisons.')
end
warning on MATLAB:polyfit:PolyNotUnique

% --------------------------------------------------------------------
function LDMatrixPlot_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LDDistancePlot_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function i_try2dolinkdisequmatrix(aln,txt,outn)

try
   warning off MATLAB:polyfit:PolyNotUnique
   outputmatrix=1; warninglarge=1;
   [D_raw,D_prime,R2,R,D_dist,siteidx,done] = linkdisequ(aln,outputmatrix,warninglarge);
   warning on MATLAB:polyfit:PolyNotUnique
   if ~(done), return; end
catch
    warndlg('No pairwise comparisons.')
    warning on MATLAB:polyfit:PolyNotUnique
    return;
end
switch (outn)
    case (1)
         Data=D_raw;
    case (3)
         Data=D_prime;
    case (6)
         Data=R2;
    case (2)
        Data=abs(D_raw);
    case (4)
        Data=abs(D_prime);
    case (5)
        Data=R;
end

	figure;
	imagesc(Data); axis equal;
	colormap(copper);
	colorbar; xlabel(txt); ylabel(txt);




% --------------------------------------------------------------------
function LDMatrixPlotD_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
outputmatrix=1;
i_try2dolinkdisequmatrix(aln,'D',1);

% --------------------------------------------------------------------
function LDMatrixPlotDAbs_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotDAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequmatrix(aln,'|D|',2);

% --------------------------------------------------------------------
function LDMatrixPlotDPrime_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotDPrime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequmatrix(aln,'D''',3);

% --------------------------------------------------------------------
function LDMatrixPlotDPrimeAbs_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotDPrimeAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequmatrix(aln,'|D''|',4);


% --------------------------------------------------------------------
function LDMatrixPlotR_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequmatrix(aln,'R',5);


% --------------------------------------------------------------------
function LDMatrixPlotR2_Callback(hObject, eventdata, handles)
% hObject    handle to LDMatrixPlotR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequmatrix(aln,'R^2',6);


%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function i_try2dolinkdisequdistance(aln,txt,outn)
try
   warning off MATLAB:polyfit:PolyNotUnique
   outputmatrix=0; warninglarge=1;
   [D_raw,D_prime,R2,R,D_dist,siteidx,done] = linkdisequ(aln,outputmatrix,warninglarge);
   warning on MATLAB:polyfit:PolyNotUnique
   if ~(done), return; end
catch
    warndlg('No pairwise comparisons.')
    warning on MATLAB:polyfit:PolyNotUnique
    return;
end
switch (outn)
    case (1)
         Data=D_raw;
    case (3)
         Data=D_prime;
    case (6)
         Data=R2;
    case (2)
        Data=abs(D_raw);
    case (4)
        Data=abs(D_prime);
    case (5)
        Data=R;
end

	figure;
	plot(D_dist, Data, '.')
	xlabel('Nucleotide distance')
	ylabel(txt)


% --------------------------------------------------------------------
function LDDistancePlotD_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'D',1);

% --------------------------------------------------------------------
function LDDistancePlotDPrime_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotDPrime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'D''',2);

% --------------------------------------------------------------------
function LDDistancePlotR2_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'R^2',6);


% --------------------------------------------------------------------
function LDDistancePlotDAbs_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotDAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'|D|',2);



% --------------------------------------------------------------------
function LDDistancePlotDPrimeAbs_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotDPrimeAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'|D''|',4);


% --------------------------------------------------------------------
function LDDistancePlotR_Callback(hObject, eventdata, handles)
% hObject    handle to LDDistancePlotR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
i_try2dolinkdisequdistance(aln,'R',5);


% --------------------------------------------------------------------
function IncludeExcludeSites_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global aln aln_ori
   prompt={'From site:','to:'};
   m=length(aln.seq);
   def={'1',num2str(m)};
   dlgTitle='Region to Analyse';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
      nstart=str2num(answer{1});
      nend=str2num(answer{2});
   if ~(nstart==1 && nend==m),
    aln_ori=aln;
    aln.seq=aln.seq(:,[nstart:nend]);
    aln.pos=aln.pos([nstart:nend]);
    viewseq(aln);
    SetMenuStatus(handles);
   end
end


% --------------------------------------------------------------------
function TestOfIndependence2x2Table_Callback(hObject, eventdata, handles)
% hObject    handle to TestOfIndependence2x2Table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    ds=189;
    ps=113;
    dr=42;
    pr=10;

    gcbox=dprstable;
	handles = guihandles(gcbox);
	set(handles.ds, 'String', num2str(ds));
	set(handles.ps, 'String', num2str(ps));
	set(handles.dr, 'String', num2str(dr));
	set(handles.pr, 'String', num2str(pr));
	set(handles.dsps,'String',num2str(ds+ps));
	set(handles.drpr,'String',num2str(dr+pr));
	set(handles.dsdr,'String',num2str(ds+dr));
	set(handles.pspr,'String',num2str(ps+pr));
	set(handles.total,'String',num2str(ps+pr+ds+dr));


% --------------------------------------------------------------------
function SNPTool_Callback(hObject, eventdata, handles)
% hObject    handle to SNPTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

snptool


% --------------------------------------------------------------------
function AssignPopulation_Callback(hObject, eventdata, handles)
% hObject    handle to AssignPopulations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[taxset2,v] = choosebox('Name','Interspecies Sequences','PromptString',...
'Sequences from species 1:','SelectString','Sequences from species 2:',...
'ListString',aln.seqnames);
if (v~=1),
	return;
else
    aln.population(taxset2)=max(aln.population)+1;
end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function AssignLocus_Callback(hObject, eventdata, handles)
% hObject    handle to AssignLoci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[locus2,v] = choosebox('Name','Assign Locus','PromptString',...
'Sequences from locus 1:','SelectString','Sequences from locus 2:',...
'ListString',aln.seqnames);

if (v~=1),
	return;
else
    aln.locus(locus2)=max(aln.locus)+1;
end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function RatioOfAdaptiveNonsynSub_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MKRPFTest_Callback(hObject, eventdata, handles)
% hObject    handle to MKRPFTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MKRPFTestExampleInput_Callback(hObject, eventdata, handles)
% hObject    handle to MKRPFTestExampleInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MKRPFTestWriteInputLine_Callback(hObject, eventdata, handles)
% hObject    handle to MKRPFTestWriteInputLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RatioOfAdaptiveNonsynSubSEW_Callback(hObject, eventdata, handles)
% hObject    handle to RatioOfAdaptiveNonsynSubSEW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D=sewfww_sampledata;
if isempty(D)
    return;
else
try
    [sew] = sewfww(D);
    i_dispheader('alpha, Smith and Eyre-Walker (2002)')
    fprintf('alpha=%.3f\n',sew);
    i_dispfooter
catch ME
    errordlg(ME.message)
end
end


% --------------------------------------------------------------------
function RatioOfAdaptiveNonsynSubFWW_Callback(hObject, eventdata, handles)
% hObject    handle to RatioOfAdaptiveNonsynSubFWW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D=sewfww_sampledata;
if isempty(D)
    return;
else
    try
    [~,fww] = sewfww(D);
    i_dispheader('alpha, Fay, Wycoff and Wu (2001)')
    fprintf('alpha=%.3f\n',fww);
    i_dispfooter
catch ME
    errordlg(ME.message)
end
end


% --------------------------------------------------------------------
function UnfoldedSiteFrequencySpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to SiteFrequencySpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
      str=aln.seqnames;
      [s,v] = listdlg('PromptString','Select an ancestral reference:',...
                      'SelectionMode','single',...
                      'ListString',str);
if v>0
    fprintf('NOTE: Sequence %d, ''%s'', was used as ancestral\nreference to count external-mutations.\n',...
            s,deblank(str{s}));
    try
    sfs(aln,s);
    f=sfs(aln,s);
    figure;
    bar(f)
catch ME
    errordlg(ME.message)
end
end



% --------------------------------------------------------------------
function FayWycoffWu01_Callback(hObject, eventdata, handles)
% hObject    handle to FayWycoffWu01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[~,data]=line4sewfww(aln);
D=sewfww_sampledata(data);
if isempty(D)
    return;
else
    [~,fww] = sewfww(D);
    i_dispheader('alpha, Fay, Wycoff and Wu (2001)')
    fprintf('alpha=%.3f\n',fww);
    i_dispfooter
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function SmithEyreWalker02_Callback(hObject, eventdata, handles)
% hObject    handle to SmithEyreWalker02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[~,data]=line4sewfww(aln);
D=sewfww_sampledata(data);
if isempty(D)
    return;
else
    sew=sewfww(D);
    i_dispheader('alpha, Smith and Eyre-Walker (2002)')
    fprintf('alpha=%.3f\n',sew);
    i_dispfooter
end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function FayWusTestH_Callback(hObject, eventdata, handles)
% hObject    handle to FayWusTestH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
      str=aln.seqnames;
      [s,v] = listdlg('PromptString','Select an ancestral reference:',...
                      'SelectionMode','single',...
                      'ListString',str);
if v>0
    fprintf('NOTE: Sequence %d, ''%s'', was used as ancestral\nreference to count external-mutations.\n',...
        s,deblank(str{s}));
    try
        seq=aln.seq;
        n=size(seq,1);
        idx=setdiff([1:n],s);
        ancseq=seq(s,:);
        seq=seq(idx,:);
        faywu00h_test(seq,ancseq);

    catch ME
    errordlg(ME)
    end
end

%      str=aln.seqnames;
%      [s,v] = listdlg('PromptString','Select an ancestral reference:',...
%                      'SelectionMode','single',...
%                      'ListString',str);
%if v>0
%    fprintf('NOTE: Sequence %d, ''%s'', was used as ancestral\nreference to count external-mutations.',s,deblank(str{s})))
%try
%    faywu00h_test(aln,s);
%catch
%    errordlg(lasterr)
%end
%end


% --------------------------------------------------------------------
function HaplotypeDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    hapdiv_test(aln.seq);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function FoldedSiteFrequencySpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to FoldedSiteFrequencySpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    sfs(aln);
    f=sfs(aln);
    figure;
    bar(f)
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function SiteFreqSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to SiteFreqSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Recombination_Callback(hObject, eventdata, handles)
% hObject    handle to Recombination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FourGameteTest_Callback(hObject, eventdata, handles)
% hObject    handle to MinimumNumOfRecombinationEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    seq=aln.seq;
    hudsonkaplan85rm(seq);
catch ME
    errordlg(ME.message)
end
