function pge_viewdata(aln)
%PGE_VIEWDATA - view sequence data and information

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

    ColNames = {'Name' 'Locus' 'Population' 'Count' 'Seq'};
    ColWidth = [100 60 60 60 200];

    seq=aln.seq;
    [n,m]=size(seq);

    if ~isfield(aln,'locus')
            aln.locus=ones(n,1);
    end

    if ~isfield(aln,'population')
    	aln.population=ones(n,1);

    end
    if ~isfield(aln,'count')
	aln.count=ones(n,1);
    end


    zz=cell(n,5);
    zz(:,1) = aln.seqnames;
    zz(:,2) = num2cell(aln.locus);
    zz(:,3) = num2cell(aln.population);
    zz(:,4) = num2cell(aln.count);

    X='ACGT-';
    s=X(seq);
    for (k=1:n),
	 if (m<=25)
         scell{k}=[s(k,[1:min(25,m)])];
	 else
         scell{k}=[s(k,[1:min(25,m)]),'...'];
	 end
    end
    zz(:,5) = scell;
%    tableGUI('FigName','Your Data','array',zz,'ColNames',ColNames,'ColWidth',ColWidth,'NumRows',min(10,n),...
%        'RowNumbers','y');

f = figure; % ('Position',[100 100 300 150]);
dat=zz;
columnname = ColNames; % {'Rate', 'Amount', 'Available', 'Fixed/Adj'};
columnformat = {'char','numeric','numeric','numeric','char'};
columneditable = [false false false false false];
t = uitable('Units','normalized','Position',...
            [0.1 0.1 0.9 0.9], 'Data', dat,...
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable,...
            'ColumnWidth',{'auto' 50 50 50 200});
