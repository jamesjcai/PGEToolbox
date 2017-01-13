function varargout = gunzip_nojava(files,varargin)
%GUNZIP Uncompress GNU zip files.
%
%   GUNZIP(FILES) uncompresses GNU zip files from the list of files
%   specified in FILES.  Directories recursively gunzip all of their
%   content.  The output gunzipped files have the same name, excluding the
%   extension '.gz', and are written to the same directory as the input
%   files.
%
%   FILES is a string or cell array of strings that specify the files or
%   directories to be uncompressed.  Individual files that are on the
%   MATLABPATH can be specified as partial pathnames. Otherwise an
%   individual file can be specified relative to the current directory or
%   with an absolute path. Directories must be specified relative to the
%   current directory or with absolute paths.  On UNIX systems, directories
%   may also start with a "~/" or a "~username/", which expands to the
%   current user's home directory or the specified user's home directory,
%   respectively.  The wildcard character '*' may be used when specifying
%   files or directories, except when relying on the MATLABPATH to resolve
%   a filename or partial pathname.
%
%   GUNZIP(FILES, OUTPUTDIR) writes the gunzipped file into the directory
%   OUTPUTDIR. OUTPUTDIR is created if it does not exist.
%
%   GUNZIP(URL, ...) extracts the gzip contents from an Internet URL. The
%   URL must include the protocol type (e.g., "http://"). The URL is
%   downloaded to the temp directory and deleted.
%
%   FILENAMES = GUNZIP(...) gunzips the files and returns the relative path
%   names of the gunzipped files into the string cell array, FILENAMES.
%
%   Examples
%   --------
%   % gunzip all *.gz files in the current directory
%   gunzip('*.gz');
%
%   % gunzip Cleve Moler's Numerical Computing with MATLAB examples
%   % to the output directory 'ncm'.
%   url ='http://www.mathworks.com/moler/ncm.tar.gz';
%   gunzip(url,'ncm')
%   untar('ncm/ncm.tar','ncm')
%
%   See also GZIP, TAR, UNTAR, UNZIP, ZIP.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


narginchk(1,2);
nargoutchk(0,1);

% rootDir is always ''
dirs = {'',varargin{:}};

% Check input arguments.
[files, rootDir, outputDir,dirCreated] = checkFilesDirInputs(mfilename, files, dirs{:});

try
% Check files input for URL .
[files, url, urlFilename] = checkFilesURLInput(files, {'gz'},'FILES',mfilename);

if ~url
   % Get and gunzip the files
   entries = getArchiveEntries(files, rootDir, mfilename);
   names = gunzipEntries(entries, outputDir);
else
   % Gunzip the URL
   names = gunzipURL(files{1}, outputDir, urlFilename);
end

% Return the names if requested
if nargout == 1
   varargout{1} = names;
end
catch exception
    if ~isempty(dirCreated)
        rmdir(dirCreated, 's');
    end
    rethrow(exception);
end
%--------------------------------------------------------------------------
function [files, url, urlFilename] = ...
   checkFilesURLInput(inputFiles, validExtensions, argName, fcnName)

% Assign the default return values
files = inputFiles;  
url = false;
urlFilename = '';

if numel(inputFiles) == 1 && isempty(strfind(inputFiles{1}, '*')) && ...
    ~isdir(inputFiles{1})

   % Check for a URL in the filename and for the file's existence
   [fullFileName, url] = checkfilename(inputFiles{1}, validExtensions, fcnName, ...
                                argName,true, tempdir);
   if url
     % Remove extension
     [~, urlFilename, ext] = fileparts(inputFiles{1});
     if ~any(strcmp(ext,{'.tgz','.gz'}))
        % Add the extension if the URL file is not .gz or .tgz
        % The URL may not be a GZIPPED file, but let pass
        urlFilename = [urlFilename ext];
     end
     files = {fullFileName};
   end
end


%--------------------------------------------------------------------------
function names = gunzipEntries(entries, outputDir)
streamCopier = ...
   com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
names = cell(1,numel(entries));
for k=1:numel(entries)
  [path, baseName] = fileparts(entries(k).file);
  
  % Set the outputDir to the file's path if outputDir is empty
  % and the path is not the current directory, (since relative
  % paths are returned in names).
  if isempty(outputDir) && ~isequal(strrep(path,'/',filesep), pwd) 
    outputDir = path;
  end
  names{k} = gunzipwrite(entries(k).file, outputDir, baseName, streamCopier);
end

%--------------------------------------------------------------------------
function names = gunzipURL(filename, outputDir, urlFilename)

streamCopier = ...
   com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
try
   names{1} = gunzipwrite(filename, outputDir, urlFilename, streamCopier);
   % Filename is temporary for URL
   delete(filename);
catch exception
   if ~isequal('MATLAB:gunzip:notGzipFormat', exception.identifier)
      delete(filename);
      throw(exception);
   else
      names{1} = fullfile(outputDir, urlFilename);
      if exist(names{1},'file') == 2
         delete(filename);
         error(message('MATLAB:gunzip:urlFileExists', names{ 1 }));
      else
         copyfile(filename, names{1})
         delete(filename);
      end
   end
end

%--------------------------------------------------------------------------
function gunzipFilename = gunzipwrite(gzipFilename, outputDir, baseName, streamCopier)
% GUNZIPWRITE Write a file in GNU zip format.
%
%   GUNZIPWRITE writes the file GZIPFILENAME in GNU zip format. 
%   OUTPUTDIR is the name of the directory for the output file. 
%   BASENAME is the base name of the output file.
%   STREAMCOPIER is a Java copy stream object. 
%
%   The output GUNZIPFILENAME is the full filename of the GNU unzipped file.

% Create the output filename from [outputDir baseName]
gunzipFilename = fullfile(outputDir,baseName);

% Create Java input stream from the gzipped filename.
fileInStream = [];
try
   fileInStream = java.io.FileInputStream(java.io.File(gzipFilename));
catch exception
   % Unable to access the gzipped file.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   error(message('MATLAB:gunzip:javaOpenError', gzipFilename));
end

% Create a Java GZPIP input stream from the file input stream.
try
   gzipInStream = java.util.zip.GZIPInputStream( fileInStream );
catch exception
   % The file is not in gzip format.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   error(message('MATLAB:gunzip:notGzipFormat', gzipFilename));
end

% Create a Java output stream from the input GZIP stream.
outStream = [];
try
   javaFile  = java.io.File(gunzipFilename);
   outStream = java.io.FileOutputStream(javaFile);
catch exception
   cleanup(gunzipFilename, outStream, gzipInStream, fileInStream);
   error(message('MATLAB:gunzip:javaOutputOpenError', gunzipFilename));
end

% Gunzip the file using the streamCopier.
try   
   streamCopier.copyStream(gzipInStream,outStream);
catch exception
   cleanup(gunzipFilename, outStream, gzipInStream, fileInStream);   
   error(message('MATLAB:gunzip:javaCopyStreamError', gunzipFilename));
end

% Cleanup and close the streams.
outStream.close;
gzipInStream.close;
fileInStream.close;

%--------------------------------------------------------------------------
function cleanup(filename, varargin)
% Close the Java streams in varargin and delete the filename.

% Close the Java streams.
for k=1:numel(varargin)
   if ~isempty(varargin{k})
      varargin{k}.close;
   end
end

% Delete the filename if it exists.
w = warning;
warning('off','MATLAB:DELETE:FileNotFound');
delete(filename);
warning(w);
