%{
function [genodata,markinfo] = snp_readmatrix(filename,nhead,ncol)
%SNP_READMATRIX - Reads genotype data in text matrix format

if nargin < 3, ncol=7; end
if nargin < 2, nhead=1; end

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.txt', 'Matrix Files (*.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Matrix file');
	if ~(filename), genodata=[]; markinfo=[]; return; end
	filename=[pathname,filename];
end
%}
%%

nhead=1;
ncol=7;
filename='C:\\Documents and Settings\\jamescai\\Desktop\\Goode\\All_Individ_Genotypes_DAFonly.txt';
[a1,a2,a3,a4,a5,a6,a7,marker] = textread(filename,...
	[repmat('%s',1,ncol),'%[ a-zA-Z\t]%*[^\n]'],...
	'delimiter','\t','commentstyle','shell','headerlines',nhead);



