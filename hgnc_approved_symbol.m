function [s]=hgnc_approved_symbol(uselocal)
%HGNC_APPROVED_SYMBOL - HGNC Custom Downloads
%
% [s]=hgnc_approved_symbol(uselocal)
%
%SEE ALSO:

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

if nargin<1&&exist('hgnc_approved_symbol.mat','file')
    uselocal=true;
end

if ~uselocal
    [t,status]=urlread('http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_app_sym&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit');
    if status==1
        s=textscan(t,'%s','delimiter','\n');        
        s=s{1}(2:end);        
    end
else
    load('hgnc_approved_symbol','s');
    return;
end
