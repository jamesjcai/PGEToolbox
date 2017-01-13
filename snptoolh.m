function varargout = snptoolh(varargin)
% snptoolh M-file for snptoolh.fig
%      snptoolh, by itself, creates a new snptoolh or raises the existing
%      singleton*.
%
%      H = snptoolh returns the handle to a new snptoolh or the handle to
%      the existing singleton*.
%
%      snptoolh('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in snptoolh.M with the given input arguments.
%
%      snptoolh('Property','Value',...) creates a new snptoolh or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before snptoolh_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to snptoolh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-04-07 11:23:42 -0500 (Sun, 07 Apr 2013) $
% $LastChangedRevision: 523 $
% $LastChangedBy: jcai $

% Edit the above text to modify the response to help snptoolh

% Last Modified by GUIDE v2.5 29-Dec-2012 18:29:56
snpversionstr='1.1';

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @snptoolh_OpeningFcn, ...
                   'gui_OutputFcn',  @snptoolh_OutputFcn, ...
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


% --- Executes just before snptoolh is made visible.
function snptoolh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to snptoolh (see VARARGIN)

% Choose default command line output for snptoolh
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

% UIWAIT makes snptoolh wait for user response (see UIRESUME)
% uiwait(handles.SNPTool);

%cdsnp;
%disp(sprintf('Current working directory is: %s', pwd))
%cdpge;
SetMenuStatus(handles);


% --- Outputs from this function are returned to the command line.
function varargout = snptoolh_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function SNPTool_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNPTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%snptoolh('SNPTool_CreateFcn',gcbo,[],guidata(gcbo))
movegui(hObject)

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CloseDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to CloseDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global haplodata hmarkinfo
SetMenuStatus(handles);

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
clear haplodata hmarkinfo
%cdpge;
close;

% --------------------------------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
info{1}='SNP Tool';
helpdlg(info,'About');



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%
function DisplayMarkerInfo(p,rsid,title)
	i_dispheader(title)
	%fprintf(['%5s\t%5s\n'],'Name','Value');
	for k=1:length(rsid)
	    fprintf('  %s\t%f\n',char(rsid(k)),p(k));
	end
	i_dispfooter




%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [filename]=i_savefile(s0,txta,txtb,txtc)
    [filename, pathname,filterindex] = uiputfile( ...
	       {txta, txtb;
		'*.*',  'All Files (*.*)'}, ...
		'Save as');
		if ~(filename),
		       helpdlg('File save canceled')
		return;
		end
		filename=[pathname,filename];
		if (filterindex==1)
			if (isempty(find(filename=='.')))
			filename=[filename,txtc];
			end
		end

	       [fid,Msg] = fopen(filename,'wt');
	       if fid == -1, error(Msg); end
	       fprintf(fid,'%s',s0);
	       fclose(fid);
   	       helpdlg('File saved')

%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [filename]=savehapmapfile(s0)
	txta='*.hapmap;*.hmp';
	txtb='HapMap Format Files (*.hapmap, *.hmp)';
	txtc='.hmp';
	[filename]=i_savefile(s0,txta,txtb,txtc);



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [filename]=savemarkerfile(s0)
   	    [filename, pathname,filterindex] = uiputfile( ...
	       {'*.marker;*.mrk;*.txt', 'Marker Information Files (*.marker, *.mrk, *.txt)';
		'*.*',  'All Files (*.*)'}, ...
		'Save as');
		if ~(filename),
		       helpdlg('File save canceled')
		return;
		end
		filename=[pathname,filename];
		if (filterindex==1)
			if (isempty(find(filename=='.')))
			filename=[filename,'.mrk'];
			end
		end

	       [fid,Msg] = fopen(filename,'wt');
	       if fid == -1, error(Msg); end
	       fprintf(fid,'%s',s0);
	       fclose(fid);
   	       helpdlg('File saved')



% --------------------------------------------------------------------
function OpenPhasedHaplotypeFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenPhasedHaplotypeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
try
	[haplodata,hmarkinfo]=snp_readhaplotype;
	SetMenuStatus(handles);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end



% --------------------------------------------------------------------
function DownloadPhasedHaplotypeDataFromHapmap_Callback(hObject, eventdata, handles)
% hObject    handle to DownloadPhasedHaplotypeDataFromHapmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo

[query,popcode]=selectMarkerPopcode('Chr9:660000..760000');
if (isempty(query)||isempty(popcode)), return; end

try
[hap,hmark,filename] = snp_downloadhaplotype(query,popcode,true);
catch exception
    if strcmp('MATLAB:dataread:TroubleReading',exception.identifier)
        % disp(exception.identifier);
        errordlg('Download error! Check the Internet connection and try again.')
    end
    rethrow(exception)
end
if (isempty(hap)),
   warndlg('Please check the Internet connection and try again.','Download Error');
else
   ButtonName=questdlg('Open downloaded data?', ...
                       'Open & View File', ...
                       'Open','View File','Cancel','Open');
   switch ButtonName,
     case 'Open',
         haplodata=hap;
         hmarkinfo=hmark;
     case 'View File',
    	 edit(filename);             
   end % switch
   SetMenuStatus(handles);
end


%{
if (isempty(s0)),
      warndlg('Data could not be found at www.hapmap.org.','Download Error');
else
   ButtonName=questdlg('Save downloaded data?', ...
                       'Save File', ...
                       'Save & Open','Just Open','Just Save','Just Open');
   switch ButtonName,
     case 'Save & Open',
    	filename=savehaplotypefile(s0);
        if ~(isempty(filename)),
            try
               [haplodata,hmarkinfo]=snp_readhaplotype(filename);
            catch exception
                errordlg('Internal Error.')
            end
        end
     case 'Just Save',
	    savehaplotypefile(s0);

     case 'Just Open',
       	filename=tempname;
    	if ~(isempty(filename)),
           [fid,Msg] = fopen(filename,'wt');
	       if fid == -1, error(Msg); end
	       fprintf(fid,'%s',s0);
	       fclose(fid);
           [haplodata,hmarkinfo]=snp_readhaplotype(filename);
        end
      case 'Don''t Save',
        helpdlg('Save cancelled')
   end % switch
   SetMenuStatus(handles);
end
%}

%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [filename]=savehaplotypefile(s0)
	txta='*.phased;*.phs';
	txtb='Phased Haplotype Format Files (*.phased, *.phs)';
	txtc='.phs';
	[filename]=i_savefile(s0,txta,txtb,txtc);

% --------------------------------------------------------------------
function ExtendedHaplotypeHomozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to ExtendedHaplotypeHomozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo

%{
      nmarker=size(haplodata,2);
      n1=floor(nmarker/2);
      n2=n1;

   prompt={'Position of first core marker:','Position of last core marker:'};
   def={num2str(max(1,n1)),num2str(max(1,n2))};
   dlgTitle=sprintf('Specify the core region, 1-%d',nmarker);
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
   n1=str2num(answer{1});
   n2=str2num(answer{2});
 else
   return;
 end

 if (n1>nmarker || n2>nmarker || n1<1 || n2<1 || n2<n1)
    errordlg('Wrong core definition')
    return;
 end
%}

[n1,b] = listdlg('Name','Markers',...
                       'PromptString','Select a marker as core:',...
                       'SelectionMode','single',...
                       'ListString',hmarkinfo.rsid);
if b
    fprintf('Marker %d (%s) has been selected as core marker.\n',...
        n1,hmarkinfo.rsid{n1});
    try
	    snp_ehh(haplodata,hmarkinfo.pos,n1);
        
    catch exception
        errordlg('Internal Error.')
        rethrow(exception)
    end
end




% --------------------------------------------------------------------
function IntegratedHaplotypeScoreiHSScatterPlot_Callback(hObject, eventdata, handles)
% hObject    handle to IntegratedHaplotypeScoreiHSScatterPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo

if size(haplodata,2)>100
      answer=questdlg(sprintf('This operation may take several\nminutes. Do you want to continue?'),...
			    'Running Time Warning',...
			    'Continue','Cancel','Continue');
      switch lower(answer)
       case 'cancel'
           helpdlg('Action cancelled.')
           return;
      end
end

    try
        snp_ihsscatter(haplodata,hmarkinfo);
    catch exception
        rethrow(exception)
        errordlg('Internal Error.')
    end



% --------------------------------------------------------------------
function iHSScatterPlotExample_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
load('example_data/ihs_result_example', 'ihs', 'markinfo');
ihs=ihs;
plot(markinfo.pos,ihs,'.');
ylabel('|iHS|')
xlabel('Genomic Postion')
title('iHS Scatter Plot')
hline(2.0)




% --------------------------------------------------------------------
function SNPLocator_Callback(hObject, eventdata, handles)
% hObject    handle to SNPLocator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    snp_locator;
catch exception    
    errordlg('Internal Error.')
    rethrow(exception)    
end


% --------------------------------------------------------------------
function HapMapSNPEHH_Callback(hObject, eventdata, handles)
% hObject    handle to HapMapSNPEHH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [rsid,popcode]=selectMarkerPopcode('rs3000');
    if isempty(rsid)||isempty(popcode), return; end
    radius=500000;
    snp_ihssingle(rsid,radius,popcode,1);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end



% --------------------------------------------------------------------
function SaveAsSWEEPFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsSWEEPFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
try
    [status] =  snp_writesweep(haplodata,hmarkinfo);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function Data_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ViewHaplotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to ViewHaplotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
    NT='ACGTN';
    Seq=NT(haplodata);
	[n,slen]=size(haplodata);
	mt = 1:60:slen;
	mt = cat(1,mt',slen+1);

	for j=1:length(mt)-1
	for i=1:n
        xs=char(Seq(i,[mt(j):mt(j+1)-1]));
		if j==1
			fprintf('%10s %s\n', sprintf('haplo_%d',i), xs);
		else
			fprintf('%10s %s\n', ' ', xs);
		end
	end
		fprintf('\n');
	end


% --------------------------------------------------------------------
function VisualHaplotype_Callback(hObject, eventdata, handles)
% hObject    handle to VisualHaplotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
figure;
try
    i_dispheader('Marker Positions')
    numvlabel(hmarkinfo.pos);
    snp_vhview(haplodata);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end



% --------------------------------------------------------------------
function Utilities_Callback(hObject, eventdata, handles)
% hObject    handle to VisualHaplotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function AssignHaplotypeDataIntoWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to VisualHaplotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
errorstatus=0;

if ~isempty(haplodata),
	assignin('base','haplodata',haplodata);
else
    errorstatus=1;
    warndlg('Haplotype data is empty.','Assign varibles into workspace')
    return;
end

if ~isempty(hmarkinfo),
	assignin('base','hmarkinfo',hmarkinfo);
else
	errorstatus=1;
	warndlg('Haplotype information is empty.','Assign varibles into workspace')
	return;
end

if (errorstatus==0),
	helpdlg('HAPLODATA and HMARKINFO assigned.','Assign varibles into workspace')
end

% --------------------------------------------------------------------
function ViewHaplotypeInfo_Callback(hObject, eventdata, handles)
% hObject    handle to ViewHaplotypeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo

try
i_dispheader('View Haplotype Information')
[n,m]=size(haplodata);
fprintf('# of samples (individuals) = %d\n# of chromosomes (haplotypes) = %d\n# of markers (SNPs) = %d\n',round(n/2),n,m);
hmarkinfo
i_dispfooter
catch exception
    errordlg('Internal Error.')
    rethrow(exception)
end


% --------------------------------------------------------------------
function SaveAsRestoredSequenceFile_Callback(hObject, eventdata, handles)
% hObject    handle to ViewHaplotypeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo

try
aln.seq=haplodata;
aln.seqtype=1;
aln.seqnames=mat2cellstr(1:size(haplodata,1));
    writefasta(aln);
    helpdlg('Sequence file saved.')
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end


function SetMenuStatus(handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata

set(handles.ChromosomeSegmentHomozygosity,'enable','off');

cannotemptyhapldata = [
    handles.ExtendedHaplotypeHomozygosity,...
    handles.IntegratedHaplotypeScoreiHSScatterPlot,...
    handles.SaveAsSWEEPFile,...
    handles.SaveAsRestoredSequenceFile,...
    handles.AssignHaplotypeDataIntoWorkspace,...
    handles.CloseDataFile,...
    handles.ViewHaplotypeInfo,...
    handles.ViewHaplotypeData,...
    handles.VisualHaplotype,...
    handles.ExcludeMonomorphicSNPs,...
    handles.MAFFilter,...
    handles.IncludeExcludeMarkers,...
    handles.IncludeExcludeIndividuals,...
    handles.WorkOnGenotypeData,...
    handles.HaplotypeNucleotideDiversity,...
    handles.VisualNucleotideDiversity,...
    handles.HaplotypeEntropy,...
    handles.HaplotypeHomozygosity,...
    handles.Haplosimilarity,...
    handles.HaplotypeSharingScore,...
    handles.ChromosomeSegmentHomozygosity,...
    handles.FourGameteTest,...
    handles.LDPlot,...
    handles.ExtendedHaplotypeHomozygosity,...
    handles.IntegratedHaplotypeScoreiHSScatterPlot,...
    handles.iHSScatterPlotExample];

if (isempty(haplodata)),
    set(cannotemptyhapldata,'enable','off');
else
    set(cannotemptyhapldata,'enable','on');
end







% --------------------------------------------------------------------
function IncludeExcludeMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
[s,v] = choosebox('Name','Pick including marker(s)','PromptString',...
    'Markers available:','SelectString','Selected markers:',...
    'ListString',hmarkinfo.rsid);
if (v==1),
    haplodata = haplodata(:,s);
    hmarkinfo.rsid=hmarkinfo.rsid(s);
    %hmarkinfo.allele=hmarkinfo.allele(s);
    %hmarkinfo.strand=hmarkinfo.strand(s);
    %hmarkinfo.chr=hmarkinfo.chr(s);
    hmarkinfo.pos=hmarkinfo.pos(s);
    hmarkinfo.maf=hmarkinfo.maf(s);
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function IncludeExcludeIndividuals_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeIndividuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
   [smpln,markn]=size(haplodata);
   names={};
   for (k=1:smpln),
       names{k}=['Idv_',num2str(k)];
   end

[s,v] = choosebox('Name','Pick including individual(s)','PromptString',...
    'Individuals available:','SelectString','Selected individuals:',...
    'ListString',names);
if (v==1),
    haplodata = haplodata(s,:);
    SetMenuStatus(handles);
else
    return;
end


% --------------------------------------------------------------------
function ChromosomeSegmentHomozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to ChromosomeSegmentHomozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Under development.')


% --------------------------------------------------------------------
function ExcludeMonomorphicSNPs_Callback(hObject, eventdata, handles)
% hObject    handle to ExcludeMonomorphicSNPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global haplodata hmarkinfo
try
    pmaf=hmarkinfo.maf;
    idx=find(pmaf>0);
    haplodata=haplodata(:,idx);
    hmarkinfo.maf=hmarkinfo.maf(idx);
    hmarkinfo.pos=hmarkinfo.pos(idx);
    hmarkinfo.rsid=hmarkinfo.rsid(idx);
    disp('Monomorphic SNPs (MAF=0) have been excluded.')
    if isempty(haplodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end

% --------------------------------------------------------------------
function MAFFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MAFFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
   prompt={'Minor allele frequency (MAF [>=], 0 to 0.5)'};
   def={'0.05'};
   dlgTitle='MAF cutoff';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
try
    cutoff=str2num(answer{1});
    if ~(cutoff>=0 && cutoff<=0.5)
        error('Incorrect MAF cutoff.')
    end
    pmaf=hmarkinfo.maf;

    idx=find(pmaf>=cutoff);
    haplodata=haplodata(:,idx);
    hmarkinfo.maf=hmarkinfo.maf(idx);
    hmarkinfo.pos=hmarkinfo.pos(idx);
    hmarkinfo.rsid=hmarkinfo.rsid(idx);
    disp(sprintf('SNPs with MAF<%f have been excluded.',cutoff))
    if isempty(haplodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch exception    
    errordlg('Internal Error.')
    rethrow(exception)    
end
end




% --------------------------------------------------------------------
function FGTDMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to FourGameteTestDMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
try
    i_dispheader('Matrix of pairwise FGT test results')
    fprintf('d=0 (white), if there are less then 4 gametes\n');
    fprintf('d=1 (black), if there are 4 gametes\n');
    i_dispfooter

    h = waitbar(0.33,'Please wait...');
    [~,D]=snp_fgt(haplodata,1);
    waitbar(0.99,h);
    close(h);
    D=D+tril(-1*ones(length(D)));

    figure;
    imagesc(D);
    axis square
    %hline([0.5:length(D)],'g-')
    %vline([0.5:length(D)],'g-')
    colormap('default');
    colorindx=gray;
    colormap(colorindx([51 64,1],:));    %[53 light gray, 1 black]
    %colorbar

%    Plot for the four-gamete test. A blackened square indicates
%that the given site pair had all four possible gametic phases.
%The presence of all four gametes in site pairs also suggests intragenic
%recombination. The diagonal line, with exons 4–9 labeled, indicates
%the location of each varying site along the gene.

catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end

% --------------------------------------------------------------------
function FourGameteTest_Callback(hObject, eventdata, handles)
% hObject    handle to FourGameteTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FGTMinNumRecombinationEvents_Callback(hObject, eventdata, handles)
% hObject    handle to FourGameteMinNumRecombinationEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global haplodata
try
    snp_fgt(haplodata);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end

% --------------------------------------------------------------------
function OpenPHASEOuptputFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenPHASEOuptputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
try
	[haplodata,hmarkinfo]=snp_readphaseout;
	SetMenuStatus(handles);
catch exception
    errordlg('Internal Error.')
    rethrow(exception)    
end

% --------------------------------------------------------------------
function WorkOnGenotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnGenotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo genodata gmarkinfo
     aswer=questdlg('HAPLOTYPE data will be converted into GENOTYPE data. Do you want to continue?','HAPLODATA => GENODATA');
     switch aswer,
     case 'Yes',
        genodata=snp_hap2geno(haplodata);
        gmarkinfo=hmarkinfo;
        try
            close;
            disp('Working on genotype data (GENODATA).')
            snptoolg
        catch exception
            rethrow(exception)
            errordlg('Internal Error.')
        end
     case 'No',
     case 'Cancel',
     otherwise

     end


% --------------------------------------------------------------------
function LDPlot_Callback(hObject, eventdata, handles)
% hObject    handle to LDPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
    try
    z=figure;
    [ldinfo]=snp_ldpairh(haplodata);
    snp_ldplot(ldinfo);
    catch exception
        if z, close(z); end
        errordlg('Internal Error.')
        rethrow(exception)        
    end


% --------------------------------------------------------------------
function HaplotypeHomozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeHomozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
i_dispheader('Haplotype Homozygosity')
fprintf ('%f\n', snp_haphom(haplodata));
i_dispfooter

% --------------------------------------------------------------------
function HaplotypeSharingScore_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeSharingScore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Haplosimilarity_Callback(hObject, eventdata, handles)
% hObject    handle to Haplosimilarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
snp_hapsimilarity(haplodata,[],true);
%i_dispheader('Haplotype Similarity')
%fprintf ('%f\n',res);
%i_dispfooter

% --------------------------------------------------------------------
function HaplotypeEntropy_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeEntropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
i_dispheader('Haplotype Entropy')
fprintf ('%f\n', snp_hapentropy(haplodata));
i_dispfooter

% --------------------------------------------------------------------
function HaplotypeNucleotideDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to HaplotypeNucleotideDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
i_dispheader('Haplotype Nucleotide Diversity')
fprintf ('%f\n', snp_hapfreqxdist(haplodata));
i_dispfooter

% --------------------------------------------------------------------
function VisualNucleotideDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to VisualNucleotideDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata
figure;
snp_hapnucdivshow(haplodata);


% --------------------------------------------------------------------
function IntegratedHaplotypeScoreiHS_Callback(hObject, eventdata, handles)
% hObject    handle to IntegratedHaplotypeScoreiHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global haplodata hmarkinfo
%{
nmarker=size(haplodata,2);
      k=floor(nmarker/2);

   prompt={'Position of core marker:'};
   def={num2str(max(1,k))};
   dlgTitle=sprintf('Specify the core marker, 1-%d',nmarker);
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
   k=str2num(answer{1});
 else
   return;
 end
%}

[n1,b] = listdlg('Name','Markers',...
                       'PromptString','Select a marker as core:',...
                       'SelectionMode','single',...
                       'ListString',hmarkinfo.rsid);
if b
    try
	    snp_ihs(haplodata,hmarkinfo,n1);
    catch exception
    errordlg('SNP_IHS: Internal Error.')        
    rethrow(exception)        
	    %warndlg('Inappropriate WINSIZE or STEP.')
   end
end
