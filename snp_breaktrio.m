function [geno] = snp_breaktrio(geno, markerpopid)
%SNP_BREAKTRIO - breaks trio by removing child individuals
%
% Usage: [geno]=snp_breaktrio(geno,markerpopid)
%
%For the CEU and YRI the HapMap genotypes are arranged as follows:
%row 1 - trio 1 parent 1
%row 2 - trio 1 parent 2
%row 3 - trio 1 child
%row 4 - trio 2 parent 1
%row 5 - trio 2 parent 2
%row 6 - trio 2 child
%
%For the JPT+CHB the genotypes are arranged as
%row 1 - individual 1
%row 2 - individual 2
%row 3 - individual 3
%row 4 - individual 4
%
%See also: SNP_HAPMAPCHILD

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

n = size(geno, 1);
if n < 3, error('Genotype data is not in trio'); end

if nargin < 2
    idx = i_defaultchildidx(geno);
elseif isempty(markerpopid)
    idx = i_defaultchildidx(geno);
else
    if isstruct(markerpopid)

        if isfield(markerpopid, 'popid')
            popid = markerpopid.popid;
            switch upper(popid)
                case {'CEU', 'YRI'}
                    idx = snp_hapmapchild(popid);
                otherwise
                    idx = i_defaultchildidx(geno);
            end
        else
            [id] = selectTrioBreakMethod;
            switch id
                case 1
                    idx = snp_hapmapchild('YRI');
                case 2
                    idx = snp_hapmapchild('CEU');
                case 3
                    idx = i_defaultchildidx(geno);
                otherwise
                    error('POPID is unknown.')
            end
        end

    elseif ischar(markerpopid) % markerpopid is popid
        %       idx=snp_hapmapchild(markerpopid);
        popid = markerpopid;
        switch upper(popid)
            case {'CEU', 'YRI'}
                idx = snp_hapmapchild(popid);
            otherwise
                idx = i_defaultchildidx(geno);
        end

    else
        error('POPID is unknown.')
    end
end

if ~isempty(idx)
    if max(idx) > size(geno, 1)
        error('SNP_BREAKTRIO found no child indv.')
    else
        geno(idx, :) = [];
    end
end


    function idx = i_defaultchildidx(geno)
        n = size(geno, 1);
        idx = 3:3:n;
        %    warning('SNP_BREAKTRIO assumes the third individuals are children.')
