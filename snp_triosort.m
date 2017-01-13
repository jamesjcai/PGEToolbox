function [idx,geno]=snp_triosort(popid,geno)
%SNP_TRIOSORT - sort HapMap genotype by order of family trio (F-M-C)
%
% [idx,sortedgeno]=snp_triosort(popid,unsortedgeno)
%
% [idx]=snp_triosort('CEU')
% [idx]=snp_triosort('YRI')
% [idx,genoordered]=snp_orderbytrio('CEU',geno)
% [idx,genoordered]=snp_orderbytrio('YRI',geno)
%
% Information source: http://www.hapmap.org/downloads/samples_individuals/

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

inputYRI={'NA18500','NA18501','NA18502','NA18503','NA18504','NA18505','NA18506','NA18507','NA18508','NA18515','NA18516','NA18517','NA18521','NA18522','NA18523','NA18852','NA18853','NA18854','NA18855','NA18856','NA18857','NA18858','NA18859','NA18860','NA18861','NA18862','NA18863','NA18870','NA18871','NA18872','NA18912','NA18913','NA18914','NA19092','NA19093','NA19094','NA19098','NA19099','NA19100','NA19101','NA19102','NA19103','NA19116','NA19119','NA19120','NA19127','NA19128','NA19129','NA19130','NA19131','NA19132','NA19137','NA19138','NA19139','NA19140','NA19141','NA19142','NA19143','NA19144','NA19145','NA19152','NA19153','NA19154','NA19159','NA19160','NA19161','NA19171','NA19172','NA19173','NA19192','NA19193','NA19194','NA19200','NA19201','NA19202','NA19203','NA19204','NA19205','NA19206','NA19207','NA19208','NA19209','NA19210','NA19211','NA19221','NA19222','NA19223','NA19238','NA19239','NA19240'};
inputCEU={'NA06985','NA06991','NA06993','NA06994','NA07000','NA07019','NA07022','NA07029','NA07034','NA07048','NA07055','NA07056','NA07345','NA07348','NA07357','NA10830','NA10831','NA10835','NA10838','NA10839','NA10846','NA10847','NA10851','NA10854','NA10855','NA10856','NA10857','NA10859','NA10860','NA10861','NA10863','NA11829','NA11830','NA11831','NA11832','NA11839','NA11840','NA11881','NA11882','NA11992','NA11993','NA11994','NA11995','NA12003','NA12004','NA12005','NA12006','NA12043','NA12044','NA12056','NA12057','NA12144','NA12145','NA12146','NA12154','NA12155','NA12156','NA12234','NA12236','NA12239','NA12248','NA12249','NA12264','NA12707','NA12716','NA12717','NA12740','NA12750','NA12751','NA12752','NA12753','NA12760','NA12761','NA12762','NA12763','NA12801','NA12802','NA12812','NA12813','NA12814','NA12815','NA12864','NA12865','NA12872','NA12873','NA12874','NA12875','NA12878','NA12891','NA12892'};


% downloaed form
% outputYRI={'NA18500','NA18502','NA18501','NA18503','NA18505','NA18504','NA18506','NA18508','NA18507','NA18860','NA18858','NA18859','NA18515','NA18517','NA18516','NA18521','NA18523','NA18522','NA18872','NA18870','NA18871','NA18854','NA18852','NA18853','NA18857','NA18855','NA18856','NA18863','NA18861','NA18862','NA18914','NA18912','NA18913','NA19094','NA19093','NA19092','NA19103','NA19102','NA19101','NA19139','NA19137','NA19138','NA19202','NA19201','NA19200','NA19173','NA19172','NA19171','NA19205','NA19204','NA19203','NA19211','NA19209','NA19210','NA19208','NA19206','NA19207','NA19161','NA19159','NA19160','NA19221','NA19222','NA19223','NA19120','NA19116','NA19119','NA19142','NA19140','NA19141','NA19154','NA19152','NA19153','NA19145','NA19143','NA19144','NA19129','NA19127','NA19128','NA19132','NA19131','NA19130','NA19100','NA19099','NA19098','NA19194','NA19193','NA19192','NA19240','NA19238','NA19239'};
% disp('YRI trios ordered as: child-mother-father');
outputYRI={'NA18501','NA18502','NA18500','NA18504','NA18505','NA18503','NA18507','NA18508','NA18506','NA18859','NA18858','NA18860','NA18516','NA18517','NA18515','NA18522','NA18523','NA18521','NA18871','NA18870','NA18872','NA18853','NA18852','NA18854','NA18856','NA18855','NA18857','NA18862','NA18861','NA18863','NA18913','NA18912','NA18914','NA19092','NA19093','NA19094','NA19101','NA19102','NA19103','NA19138','NA19137','NA19139','NA19200','NA19201','NA19202','NA19171','NA19172','NA19173','NA19203','NA19204','NA19205','NA19210','NA19209','NA19211','NA19207','NA19206','NA19208','NA19160','NA19159','NA19161','NA19223','NA19222','NA19221','NA19119','NA19116','NA19120','NA19141','NA19140','NA19142','NA19153','NA19152','NA19154','NA19144','NA19143','NA19145','NA19128','NA19127','NA19129','NA19130','NA19131','NA19132','NA19098','NA19099','NA19100','NA19192','NA19193','NA19194','NA19239','NA19238','NA19240'};
outputCEU={'NA12003','NA12004','NA10838','NA12005','NA12006','NA10839','NA12056','NA12057','NA10851','NA11829','NA11830','NA10856','NA11831','NA11832','NA10855','NA11992','NA11993','NA10860','NA11994','NA11995','NA10861','NA12154','NA12236','NA10830','NA12155','NA12156','NA10831','NA07357','NA07345','NA07348','NA06994','NA07000','NA07029','NA07022','NA07056','NA07019','NA11839','NA11840','NA10854','NA12043','NA12044','NA10857','NA12144','NA12145','NA10846','NA12146','NA12239','NA10847','NA12264','NA12234','NA10863','NA12716','NA12717','NA12707','NA12891','NA12892','NA12878','NA12812','NA12813','NA12801','NA12814','NA12815','NA12802','NA12872','NA12873','NA12864','NA12874','NA12875','NA12865','NA12760','NA12761','NA12752','NA12762','NA12763','NA12753','NA07034','NA07055','NA07048','NA06993','NA06985','NA06991','NA12750','NA12751','NA12740','NA11881','NA11882','NA10859','NA12248','NA12249','NA10835'};

switch upper(popid)
    case 'YRI'
        [y,idx]=ismember(inputYRI,outputYRI);
        disp('YRI trios ordered as: father-mother-child');
    case 'CEU'
        [y,idx]=ismember(inputCEU,outputCEU);
        disp('CEU trios ordered as: father-mother-child');
    otherwise
        error('Incorrect POPID');
end
if nargout>1
    geno=geno(idx,:);
end