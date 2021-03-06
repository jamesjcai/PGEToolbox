function varargout = snptoolg(varargin)
% snptoolg M-file for snptoolg.fig
%      snptoolg, by itself, creates a new snptoolg or raises the existing
%      singleton*.
%
%      H = snptoolg returns the handle to a new snptoolg or the handle to
%      the existing singleton*.
%
%      snptoolg('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in snptoolg.M with the given input arguments.
%
%      snptoolg('Property','Value',...) creates a new snptoolg or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before snptoolg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to snptoolg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $

% Edit the above text to modify the response to help snptoolg

% Last Modified by GUIDE v2.5 28-Apr-2013 12:05:44
snpversionstr='1.1';

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @snptoolg_OpeningFcn, ...
                   'gui_OutputFcn',  @snptoolg_OutputFcn, ...
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


% --- Executes just before snptoolg is made visible.
function snptoolg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to snptoolg (see VARARGIN)

% Choose default command line output for snptoolg
handles.output = hObject;

scrsz=get(0,'ScreenSize');
pos_act=get(gcf,'Position');
xr=scrsz(3)-pos_act(3);
xp=round(xr/2);
yr=scrsz(4)-pos_act(4);
yp=round(yr/2);
set(gcf,'position',[xp yp pos_act(3) pos_act(4)]);

% if ~ispref('PGEToolbox','lastworkingdir')
%    addpref('PGEToolbox','lastworkingdir',pwd)
% else
%    try
%     cd(getpref('PGEToolbox','lastworkingdir',pwd));
%    catch ME
%     rmpref('PGEToolbox','lastworkingdir');
%    end
% end
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes snptoolg wait for user response (see UIRESUME)
% uiwait(handles.SNPToolg);
%clear global genodata gmarkinfo
%cdsnp;
%disp(sprintf('Current working directory is: %s', pwd))
%cdpge;
SetMenuStatus(handles);


% --- Outputs from this function are returned to the command line.
function varargout = snptoolg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function SNPToolg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNPToolg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%snptoolg('SNPTool_CreateFcn',gcbo,[],guidata(gcbo))
movegui(hObject)

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenHapMapGenotypeFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenHapMapGenotypeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[genodata,gmarkinfo]=snp_readhapmap;
SetMenuStatus(handles);

% --------------------------------------------------------------------
function OpenLinkageGenotypeFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenLinkageGenotypeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[genodata,gmarkinfo]=snp_readlinkage;
SetMenuStatus(handles);

% --------------------------------------------------------------------
function CloseDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to CloseDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global genodata gmarkinfo
SetMenuStatus(handles);

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
clear genodata gmarkinfo
close;

% --------------------------------------------------------------------
function Data_Callback(hObject, eventdata, handles)
% hObject    handle to Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --------------------------------------------------------------------
function AssignGenotypeDataIntoWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to AssignGenotypeDataIntoWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
errorstatus=0;

if ~isempty(genodata),
	assignin('base','genodata',genodata);
else
    errorstatus=1;
    warndlg('Genotype data is empty.','Assign varibles into workspace')
    return;
end

if ~isempty(gmarkinfo),
	assignin('base','gmarkinfo',gmarkinfo);
else
	errorstatus=1;
	warndlg('Genotype information is empty.','Assign varibles into workspace')
	return;
end

if (errorstatus==0),
	helpdlg('GENODATA and GMARKINFO assigned.','Assign varibles into workspace')
end


% --------------------------------------------------------------------
function ViewGenotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to ViewGenotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    snp_viewgeno(genodata)
catch ME
    errordlg(ME.message)
end
%viewgeno(genotranspose(genodata),'markers')

% --------------------------------------------------------------------
function MinmumAlleleFreq_Callback(hObject, eventdata, handles)
% hObject    handle to MinmumAlleleFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    p=snp_maf(genodata);
    DisplayMarkerInfo(p,gmarkinfo.rsid,'Minor Allele Frequency')
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ObservedHeterozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to ObservedHeterozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    h=snp_obshet(genodata);
    DisplayMarkerInfo(h,gmarkinfo.rsid,'Observed Heterozygosity')
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function PredictedHeterozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to PredictedHeterozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    h=snp_predhet(genodata);
    DisplayMarkerInfo(h,gmarkinfo.rsid,'Predicted Heterozygosity')
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function GenoPercentage_Callback(hObject, eventdata, handles)
% hObject    handle to GenoPercentage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    p=snp_genoprct(genodata);
    DisplayMarkerInfo(p,gmarkinfo.rsid,'% Geno')
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function HardyWeinbergPValue_Callback(hObject, eventdata, handles)
% hObject    handle to HardyWeinbergPValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
	warning off MATLAB:divideByZero
	p=snp_hwpvalue(genodata,gmarkinfo);
	warning on MATLAB:divideByZero
	%p'
	DisplayMarkerInfo(p,gmarkinfo.rsid,'Hardy-Weinberg p-value')
catch ME
    errordlg(ME.message)
end


%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%
function DisplayMarkerInfo(p,rsid,title)
	i_dispheader(title)
	%fprintf(['%5s\t%5s\n'],'Name','Value');
	for (k=1:length(rsid)),
	    fprintf(['%d\t%f\t%s\n'],k,p(k),char(rsid(k)));
	end
	i_dispfooter


% --------------------------------------------------------------------
function CheckMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to CheckMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo

try
    ColNames = {'Name' 'Position' 'ObsHET' 'PredHET' 'HWpval' '%Geno' 'MAF'};
    ColWidth = [80 100 60 60 60 60 60];

    n=length(gmarkinfo.rsid);
    zz=cell(n,7);
    zz(:,1) = gmarkinfo.rsid;
    zz(:,2) = num2cell(gmarkinfo.pos);

    h=snp_obshet(genodata);
    zz(:,3) = num2cell(h');

    h=snp_predhet(genodata);
    zz(:,4) = num2cell(h');

    p=snp_genoprct(genodata);
    zz(:,6) = num2cell(p)';

    m=snp_maf(genodata);
    zz(:,7) = num2cell(m');

    warning off MATLAB:divideByZero
    p=snp_hwpvalue(genodata);
    warning on MATLAB:divideByZero
    zz(:,5) = num2cell(p)';

    tableGUI('FigName','Check Markers','array',zz,'ColNames',ColNames,'ColWidth',ColWidth,'NumRows',12,...
             'RowNumbers','y','checks','y');

catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function DownloadHapMapGenotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to DownloadHapMapGenotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo

[marker,popcode]=selectMarkerPopcode;
if (isempty(marker)&&isempty(popcode)),
   % helpdlg('Action cancelled.','User Input')
   return;
elseif (isempty(marker)&&~isempty(popcode)),
   warndlg('MARKER (or REGION) cannot be empty.','Input Error')
   return;
end

try
    [geno,gmark,filename] = snp_downloadhapmap(marker,popcode,1);
catch exception
    if strcmp('MATLAB:dataread:TroubleReading',exception.identifier)
        % disp(exception.identifier);
        errordlg('Download error! Check the Internet connection and try again.')
    end
    return;
end
if (isempty(geno)),
   warndlg('Please check the Internet connection and try again.','Download Error');
else
   ButtonName=questdlg('Open downloaded data?', ...
                       'Open & View File', ...
                       'Open','View File','Cancel','Open');
   switch ButtonName,
     case 'Open',
         genodata=geno;
         gmarkinfo=gmark;
     case 'View File',
    	 edit(filename);             
   end % switch
   SetMenuStatus(handles);
end

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
%%% SUB ccc %%%
%%%%%%%%%%%%

function [filename]=savehapmapfile(s0,filename)
    if nargin<2
	txta='*.hapmap;*.hmp';
	txtb='HapMap Format Files (*.hapmap, *.hmp)';
	txtc='.hmp';
	[filename]=i_savefile(s0,txta,txtb,txtc);
    else
           [fid,Msg] = fopen(filename,'wt');
	       if fid == -1, error(Msg); end
	       fprintf(fid,'%s',s0);
	       fclose(fid);
    end


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
function SaveAsSNPHAPFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsSNPHAPFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    [status] =  snp_writesnphap(genodata);
catch ME
    errordlg(ME.message)
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function IncludeExcludeMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[s,v] = choosebox('Name','Pick including marker(s)','PromptString',...
    'Markers available:','SelectString','Selected markers:',...
    'ListString',gmarkinfo.rsid');
if (v==1),
    %s2=s*2-1;
    %sx=zeros(1,length(s)*2);
    %for (k=1:length(s));
    %    sx(k*2-1)=s2(k);
    %    sx(k*2)=s2(k)+1;
    %end
    %genodata = genodata(:,sx);
    %gmarkinfo.rsid=gmarkinfo.rsid(s);
    %gmarkinfo.allele=gmarkinfo.allele(s);
    %gmarkinfo.strand=gmarkinfo.strand(s);
    %gmarkinfo.chr=gmarkinfo.chr(s);
    %gmarkinfo.pos=gmarkinfo.pos(s);

    %[genodata,gmarkinfo]=snp_subgeno(s,genodata,gmarkinfo);
    [genodata,gmarkinfo]=snp_pickmarker(genodata,gmarkinfo,s);
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function IncludeExcludeIndividuals_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
   [smpln]=size(genodata,1);
   names={};
   for k=1:smpln
       names{k}=['Idv_',num2str(k)];
   end

[s,v] = choosebox('Name','Pick including individual(s)','PromptString',...
    'Individuals available:','SelectString','Selected individuals:',...
    'ListString',names);
if (v==1),
    genodata = genodata(s,:);
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function LD_Callback(hObject, eventdata, handles)
% hObject    handle to LD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AlleleFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to AlleleFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    snp_allefreq(genodata,gmarkinfo);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function GenotypeFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to GenotypeFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    snp_genofreq(genodata,gmarkinfo);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function FrequencyPie_Callback(hObject, eventdata, handles)
% hObject    handle to FrequencyPie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    try
        z=figure;
        snp_freqpie
    catch
        close(z);
        errordlg(lasterr)
   end


% --------------------------------------------------------------------
function BioMartSNP_Callback(hObject, eventdata, handles)
% hObject    handle to BioMartSNP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%[s1,gmarkinfo]=biomartsnp;
[s1]=biomartsnp;
if (isempty(s1)),
      warndlg('No data downloaded. Please check the Internet connection and try again.','Download Error');
else
   ButtonName=questdlg('Save downloaded data?', ...
                       'Save File', ...
                       'Save & Open','Just Save','Don''t Save','Save & Open');

   switch ButtonName,
     case 'Save & Open',
	filename=savemarkerfile(s1);
	if ~(isempty(filename)),
		x=inputdlg(filename,'View Saved File',10,{s1});
	end
     case 'Just Save',
	    savemarkerfile(s1);
      case 'Don''t Save',
        helpdlg('Save cancelled')
   end % switch
   SetMenuStatus(handles);
end


% --------------------------------------------------------------------
function ViewMarkerInfo_Callback(hObject, eventdata, handles)
% hObject    handle to ViewMarkerInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gmarkinfo
try
	snp_viewmark(gmarkinfo)
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function VisualGenotype_Callback(hObject, eventdata, handles)
% hObject    handle to VisualGenotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
figure;
try
    if isfield('pos',gmarkinfo)
    i_dispheader('Marker Positions')
    numvlabel(gmarkinfo.pos);
    end
    snp_vgview(genodata);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function PlotMarkerPosition_Callback(hObject, eventdata, handles)
% hObject    handle to PlotMarkerPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
	gmarkinfo.maf=snp_maf(genodata);
	z=figure;
	snp_plotmarkpos(gmarkinfo);
catch
    close(z);
    errordlg(lasterr)
end


% --------------------------------------------------------------------
function PlotMarkerPositionOnChr_Callback(hObject, eventdata, handles)
% hObject    handle to PlotMarkerPositionOnChr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
	gmarkinfo.maf=snp_maf(genodata);
	z=figure;
	snp_plotmarkpos(gmarkinfo,1);
catch
	close(z);
    errordlg(lasterr)
end





% --------------------------------------------------------------------
function ViewMarkerPosition_Callback(hObject, eventdata, handles)
% hObject    handle to ViewMarkerPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
	gmarkinfo.maf=snp_maf(genodata);
	z=figure;
	snp_plotmarkpos(gmarkinfo,0);
catch
    close(z);
    errordlg(lasterr)
end





% --------------------------------------------------------------------
function SNPFrequencySpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to SNPFrequencySpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
	z=figure;
	snp_sfs(genodata);

catch
	close(z)
    errordlg(lasterr)
end


%x=snp_maf(genodata);
%figure;
    %bar(histc(x,0:0.05:0.5-eps),'facecolor',[.4 .2 .4])
%    bar(histc(x,0:0.05:0.5-eps))
%    set(gca,'XTickLabel',0.05:0.05:0.5)
%    title(sprintf('Folded SFS'))
%myhist(x)
%xlabel('SNP frequency (MAF)')
%ylabel('Proportion of SNPs')



% --------------------------------------------------------------------
function Tajima89Test_Callback(hObject, eventdata, handles)
% hObject    handle to Tajima89Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    snp_tajima89d(genodata);
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function NucleotideDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to NucleotideDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    i_dispheader('Estimation of theta (i.e., 4Nu)')
	fprintf (['Theta-Pi         = %f (per sequence)\n'], snp_thetapi(genodata,gmarkinfo,0));
	fprintf (['                 = %f (per site)\n'], snp_thetapi(genodata,gmarkinfo,1));
    i_dispfooter
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function FaysthetaH_Callback(hObject, eventdata, handles)
% hObject    handle to FaysthetaH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    i_dispheader('Estimation of theta (i.e., 4Nu)')
	fprintf (['Theta-H         = %f (per sequence)\n'], snp_thetah(genodata,gmarkinfo,0,[],0));
	fprintf (['                = %f (per site)\n'], snp_thetah(genodata,gmarkinfo,1,[],0));
    i_dispfooter
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ThetaW_Callback(hObject, eventdata, handles)
% hObject    handle to ThetaW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try

	i_dispheader('Estimation of theta (i.e., 4Nu)')
	fprintf (['Theta-W         = %f (per sequence)\n'], snp_thetaw(genodata,gmarkinfo,0));
	fprintf (['                = %f (per site)\n'], snp_thetaw(genodata,gmarkinfo,1));
    i_dispfooter
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function FayWuTestH_Callback(hObject, eventdata, handles)
% hObject    handle to FayWuTestH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
	snp_faywu00h(genodata);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function CompositeLikelihood_Callback(hObject, eventdata, handles)
% hObject    handle to CompositeLikelihood (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
snp_complike(genodata,gmarkinfo);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function LDPlot_Callback(hObject, eventdata, handles)
% hObject    handle to LDPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SNPWeblink_Callback(hObject, eventdata, handles)
% hObject    handle to SNPWeblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gmarkinfo
try
    if isempty(gmarkinfo)
	snp_weblink;
    else
	snp_weblink(gmarkinfo)
    end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function SNPDiversity_Callback(hObject, eventdata, handles)
% hObject    handle to SNPDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
snp_diversity(genodata);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function SaveAsPHASEFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsPHASEFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[status] =  snp_writephase(genodata,gmarkinfo);
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function SNPLocator_Callback(hObject, eventdata, handles)
% hObject    handle to SNPLocator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    snp_locator;
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ExcludeMonomorphicSNPs_Callback(hObject, eventdata, handles)
% hObject    handle to ExcludeMonomorphicSNPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
try
    pmaf=snp_maf(genodata);
    idx=find(pmaf>0);
    [genodata,gmarkinfo]=snp_pickmarker(genodata,gmarkinfo,idx);

    disp('Monomorphic SNPs (MAF=0) have been excluded.')
    if isempty(genodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function MAFFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MAFFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
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
    pmaf =  snp_maf(genodata);
    idx=find(pmaf>=cutoff);
    [genodata,gmarkinfo]=snp_pickmarker(genodata,gmarkinfo,idx);
    disp(sprintf('SNPs with MAF<%f have been excluded.',cutoff))
    if isempty(genodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end
end


% --------------------------------------------------------------------
function Utilities_Callback(hObject, eventdata, handles)
% hObject    handle to SNPDiversity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SetMenuStatus(handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
%load testgeno

cannotemptygeno = [
    handles.SaveAsFDIST2File,...
    handles.SaveAsFSTATFile,...
    handles.SaveAsGENEPOPFile,...
    handles.SaveAsGDAFile,...
    handles.SaveAsLinkageFile,...
    handles.SaveAsSNPHAPFile,...
    handles.SaveAsPHASEFile,...
    handles.SaveAsSTRUCTUREFile,...
    handles.PlotMarkerPositionOnChr,...
    handles.AssignGenotypeDataIntoWorkspace,...
    handles.CloseDataFile,...
    handles.SNPFrequencySpectrum,...
    handles.ViewGenotypeData,...
    handles.VisualGenotype,...
    handles.AlleleFrequency,...
    handles.GenotypeFrequency,...
    handles.ViewMarkerInfo,...
    handles.ViewMarkerPosition,...
    handles.CheckMarkersMenu,...
    handles.MinmumAlleleFreq,...
    handles.SNPFrequencySpectrum,...
    handles.ExcludeMonomorphicSNPs,...
    handles.MAFFilter,...
    handles.GenoPercentage,...
    handles.CallRateHistogram,...
    handles.GenoPercentageFilter,...
    handles.ObservedHeterozygosity,...
    handles.PredictedHeterozygosity,...
    handles.HardyWeinbergPValue,...
    handles.HWPValueFilter,...
    handles.ExcludeChildsInTrios,...
    handles.IncludeExcludeIndividuals,...
    handles.IncludeExcludeMarkers,...
    handles.WorkOnHaplotypeData,...
    handles.NucleotideDiversity,...
    handles.ThetaW,...
    handles.FaysthetaH,...
    handles.Tajima89Test,...
    handles.FayWuTestH,...
    handles.SNPDiversity,...
    handles.SNPHeterozygosity,...
    handles.CompositeLikelihood,...
    handles.SweepFinder,...
    handles.LDPlot];
%if (isempty(genodata)), disp('Empty GENODATA'); end
%if (isempty(gmarkinfo)), disp('Empty MARKINFO'); end

if (isempty(genodata)),
    set(cannotemptygeno,'enable','off');
else
    set(cannotemptygeno,'enable','on');
end


% --------------------------------------------------------------------
function rMHHrHHTests_Callback(hObject, eventdata, handles)
% hObject    handle to rMHHrHHTests (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
snp_mhhtest;
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function FstWeirCockerham_Callback(hObject, eventdata, handles)
% hObject    handle to FstWeirCockerham (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
       [marker,popid1,popid2]=selectMarker2Popcodes('EDAR',5,1);
       if isempty(marker)
           return;
       end

disp('Downloading GENODATA from population 1....')
       %[s] = snp_downloadhapmap(marker,popid1);
       %filename=tempname;
       %[fid,Msg] = fopen(filename,'wt');
       %if fid == -1, error(Msg); end
       %fprintf(fid,'%s',s);
       %fclose(fid);
       %snp_readhapmap(filename);
       
       [geno1,mark1] = snp_downloadhapmap(marker,popid1);
       if ismember(popid2,{'CEU','YRI'})
       [geno1]=snp_breaktrio(geno1,mark1);
       end
   
disp('Downloading GENODATA from population 2....')
       %[s] = snp_downloadhapmap(marker,popid2);
       %filename=tempname;
       %[fid,Msg] = fopen(filename,'wt');
       %if fid == -1, error(Msg); end
       %fprintf(fid,'%s',s);
       %fclose(fid);
       %snp_readhapmap(filename);
       [geno2,mark2] = snp_downloadhapmap(marker,popid2);
       if ismember(popid2,{'CEU','YRI'})
       [geno2]=snp_breaktrio(geno2,mark2);
       end


disp('Extracting common markers....')
x=intersect(mark1.pos,mark2.pos);
[a]=find(ismember(mark1.pos,x));
[b]=find(ismember(mark2.pos,x));
[geno1,mark1]=snp_pickmarker(geno1,mark1,a);
[geno2,mark2]=snp_pickmarker(geno2,mark2,b);


if (isempty(geno1) || isempty(geno2))
    warndlg('No SNP downloaded.')
    return;
else
    [fv]=snp_fst(geno1,geno2);
    poptext=sprintf('%s vs. %s',popid1,popid2);
    i_dispheader(['Fst (Weir & Cockerham 1984) for SNPs, ',poptext])
    for k=1:length(mark1.rsid)
    fprintf('%s\t%f\n',mark1.rsid{k},fv(k))
    end
    i_dispfooter
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function SNPHeterozygosity_Callback(hObject, eventdata, handles)
% hObject    handle to SNPHeterozygosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    snp_heterozygosity(genodata);
catch ME
    errordlg(ME.message)
end




% --------------------------------------------------------------------
function ExcludeChildsInTrios_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeChildsInTrios (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
uiwait(helpdlg('This command applies only to HapMap CEU and YRI samples.'));

choice = questdlg('Would you like to continue?', ...
	'Exclusion of Children from HapMap Trios', ...
	'Yes','No','Cancel','No');
% Handle response
switch choice
    case 'Yes'
        dessert = true;
    case 'No'
        dessert = false;
    case 'Cancel'
        dessert = false;
end
if dessert
    try
        n1=snp_samplen(genodata);
        [genodata]=snp_breaktrio(genodata,gmarkinfo);
        n2=snp_samplen(genodata);
        if n1~=n2
            i_dispheader('Exclusion of Child Individuals from Trios')
            fprintf('Sample size before exclusion = %d\n', n1);
            fprintf('Sample size after exclusion = %d\n', n2);
            i_dispfooter
            helpdlg(sprintf('%d individuls excluded.',n1-n2))
        else
            helpdlg('No child in trios excluded.')
        end
    catch ME
        errordlg(ME.message)
    end
end



% --------------------------------------------------------------------
function SaveAsSTRUCTUREInputFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsSTRUCTUREInputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[status] =  snp_writestructure(genodata,gmarkinfo);
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end



% --------------------------------------------------------------------
function GenoPercentageFilter_Callback(hObject, eventdata, handles)
% hObject    handle to GenoPercentageFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
   prompt={'Genotyped Percentage (GenoPrct [>=], 0 to 1)'};
   def={'0.95'};
   dlgTitle='GenoPrct cutoff';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
try
    cutoff=str2num(answer{1});
    if ~(cutoff>=0 && cutoff<=1)
        error('Incorrect GenoPrct cutoff.')
    end

    p=snp_genoprct(genodata);
    idx=find(p>=cutoff);
    [genodata,gmarkinfo]=snp_pickmarker(genodata,gmarkinfo,idx);
    disp(sprintf('SNPs with GenoPrct<%f have been excluded.',cutoff))
    if isempty(genodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end
end


% --------------------------------------------------------------------
function SaveAsFDIST2File_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsFDIST2File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    [status] =  snp_writefdist2(genodata);
catch ME
    errordlg(ME.message)
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function SaveAsFSTATFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsFSTATFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    [status] =  snp_writefstat(genodata);
catch ME
    errordlg(ME.message)
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function SaveAsGENEPOPFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsGENEPOPFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    [status] =  snp_writegenepop(genodata);
catch ME
    errordlg(ME.message)
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function SaveAsGDAFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsGDAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    [status] =  snp_writegda(genodata);
catch ME
    errordlg(ME.message)
end
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end


% --------------------------------------------------------------------
function SaveDataAs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function WorkOnHaplotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnHaplotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global haplodata hmarkinfo genodata gmarkinfo

     % aswer=questdlg('GENOTYPE data will be converted into pseudo-HAPLOTYPE data. Do you want to continue?',...
     %                'GENODATA => HAPLODATA');
     aswer=questdlg('GENOTYPE data will be phased uisng program PHASE. Do you want to continue?',...
                    'GENODATA => HAPLODATA');
     switch aswer,
         case 'Yes',
           % haplodata=snp_geno2hap(genodata);
           % hmarkinfo=gmarkinfo;
            [haplodata,hmarkinfo]=snp_phaserun(genodata,gmarkinfo);
            try
                close;
                disp('Working on haplotype data (HAPLODATA).')
                snptoolh
            catch ME
    errordlg(ME.message)
            end
         case 'No',

         case 'Cancel',

         otherwise
     end


% --------------------------------------------------------------------
function CallRateHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to CallRateHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    figure;
    hist(snp_genoprct(genodata),20);
    xlabel('Call Rate');
catch ME
    errordlg(ME.message)
end





% --------------------------------------------------------------------
function HWPValueFilter_Callback(hObject, eventdata, handles)
% hObject    handle to HWPValueFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global genodata gmarkinfo
   prompt={'HW P-value cutoff (0 to 1)'};
   def={'0.001'};
   dlgTitle='P-value cutoff';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
try
    cutoff=str2num(answer{1});
    if ~(cutoff>=0 && cutoff<=1)
        error('Incorrect P-value cutoff.')
    end
    pmaf =  snp_hwpvalue(genodata);
    idx=find(pmaf>=cutoff);
    [genodata,gmarkinfo]=snp_pickmarker(genodata,gmarkinfo,idx);
    disp(sprintf('SNPs with HW P-value<%f have been excluded.',cutoff))
    if isempty(genodata)
        warndlg('All markers have been excluded.')
    end
    SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end
end


% --------------------------------------------------------------------
function LDPlotSNPLDPAIR_Callback(hObject, eventdata, handles)
% hObject    handle to LDPlotSNPLDPAIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo

if (size(genodata,2)<=100),
    try
    z=figure;
    ldinfo=snp_ldpair(genodata);
    snp_ldplot(ldinfo);
    catch
        close(z);
        errordlg(lasterr)
    end
else

[s,v] = choosebox('Name','Pick including marker(s)','PromptString',...
    'Markers available:','SelectString','Selected markers:',...
    'ListString',gmarkinfo.rsid');
if (v==1),
    if (length(s)>500)
        warndlg('SNP_LDPAIR is too slow to processes >500 SNPs.')
        return;
    end
    s2=s*2-1;
    sx=zeros(1,length(s)*2);
    for (k=1:length(s));
        sx(k*2-1)=s2(k);
        sx(k*2)=s2(k)+1;
    end
    genodata = genodata(:,sx);
    gmarkinfo.rsid=gmarkinfo.rsid(s);
    gmarkinfo.allele=gmarkinfo.allele(s);
    gmarkinfo.strand=gmarkinfo.strand(s);
    gmarkinfo.chr=gmarkinfo.chr(s);
    gmarkinfo.pos=gmarkinfo.pos(s);
    SetMenuStatus(handles);
    try
        z=figure;
        tic
        ldinfo=snp_ldpair(genodata);
        snp_ldplot(ldinfo);
        toc
    catch exception
        close(z);
        errordlg(exception.identifier)
    end
else
    return;
end
end

% --------------------------------------------------------------------
function LDPlotEMLDJava_Callback(hObject, eventdata, handles)
% hObject    handle to LDPlotEMLDJava (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
try
    ldinfo=emldrun(genodata);
    if ~isempty(ldinfo)
        z=figure;
        snp_ldplot(ldinfo);
    end
catch exception
    if z,
        close(z);
    end
    errordlg(exception.identifier);
    rmpref('pgetoolbox','emldrun_prgmdir');
    %{
    ButtonName=questdlg('Locating program?', ...
	                    'Select Directory', ...
	                    'Yes','No','Cancel','Yes');
    if strcmp(ButtonName,'Yes')
        rmpref('pgetoolbox','sweepfinderrun_dir');
    end
    %}
end    



% --------------------------------------------------------------------
function SaveAsLinkageFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsLinkageFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[status]= snp_writelinkage(genodata,gmarkinfo);
if (status==1),
    helpdlg('File saved.')
else
    warndlg('File not saved.')
end




% --------------------------------------------------------------------
function SweepFinder_Callback(hObject, eventdata, handles)
% hObject    handle to SweepFinder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
%figure;
try
    sweepfinderrun(genodata,gmarkinfo);
catch exception
    errordlg(exception.identifier);
    rmpref('pgetoolbox','sweepfinderrun_dir');
    %{
    ButtonName=questdlg('Locating program?', ...
	                    'Select Directory', ...
	                    'Yes','No','Cancel','Yes');
    if strcmp(ButtonName,'Yes')
        rmpref('pgetoolbox','sweepfinderrun_dir');
    end
    %}
end


% --------------------------------------------------------------------
function Open1000GenomesGenotypeFile_Callback(hObject, eventdata, handles)
% hObject    handle to Open1000GenomesGenotypeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
[genodata,gmarkinfo]=snp_readvcf;
SetMenuStatus(handles);


% --------------------------------------------------------------------
function HaploRegSNPAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to HaploRegSNPAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    snp_haploreginfo;
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function Download1000GenomesGenotypeData_Callback(hObject, eventdata, handles)
% hObject    handle to Download1000GenomesGenotypeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata gmarkinfo
query=select1000GenomesRegion;
noerror=false;
if ~isempty(query)    
    %{
    aaa=sprintf('A web browser will open. After the web page is completely loaded, the genotype data will be automatically downloaded.');
    choice = questdlg(aaa, ...
	'Open web browser', ...
	'Continue','Cancel','Continue');    
    switch choice
        case 'Continue'            
            dessert = true;
        case 'Cancel'            
            dessert = false;
    end
    
    if dessert
        try
            [g,m]=snp_download1000genomes(query,'ALL');
            if ~isempty(g), genodata=g; end
            if ~isempty(m), gmarkinfo=m; end
            
        catch ME
            errordlg(ME.message);
        end
    end
    %}
    %uiwait(helpdlg('______'))
        noerror=true;
        try
            [g,m]=snp_download1000genomes(query,'ALL',true);
            if ~isempty(g)
                genodata=g;
            else
                noerror=false;
            end
            if ~isempty(m)
                gmarkinfo=m;
            else
                noerror=false;
            end
            
        catch ME
            noerror=false;
            errordlg(ME.message);            
        end
end
SetMenuStatus(handles);
if noerror
    helpdlg('Genotype data is loaded.','The 1000 Genomes Project')
else
    errordlg('Empty genotype data!','The 1000 Genomes Project')
end



% --------------------------------------------------------------------
function IncludeExclude1000GenomesIndividuals_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExclude1000GenomesIndividuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
[smpln]=size(genodata,1);
if smpln~=1092
   errordlg('Incorrect number of samples in genotype data.') 
else
   [idvlist]=smplpop_mapping_1kgenomes;
    [s,v] = choosebox('Name','Pick including individual(s)','PromptString',...
        'Individuals available:','SelectString','Selected individuals:',...
        'ListString',idvlist);
    if v==1
        genodata = genodata(s,:);
        SetMenuStatus(handles);
    else
        return;
    end
end


% --------------------------------------------------------------------
function IncludeExclude1000GenomesPopulations_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExclude1000GenomesPopulations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global genodata
[smpln]=size(genodata,1);
if smpln~=1092
   errordlg('Incorrect number of samples in genotype data.') 
else
    [~,poplist]=smplpop_mapping_1kgenomes;
    [upoplist,~,c]=unique(poplist);
    [s,v] = choosebox('Name','Pick including individual(s)','PromptString',...
        'Individuals available:','SelectString','Selected individuals:',...
        'ListString',upoplist');
    i=ismember(c,s);
    if v==1
        genodata=genodata(i,:);
        SetMenuStatus(handles);
    else
        return;
    end
end
