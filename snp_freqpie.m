function snp_freqpie(ishapmap3,marker)
%SNP_FREQPIE - Pie chart of allele and genotype frequencies among populations
%Usage: snp_freqpie('rs2693665')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<1),
    ishapmap3=0;
end

if (nargin<2),
   prompt={sprintf('This command will generate a pie chart of allele and genotype frequencies.\n Please enter a genotyped SNP ID: ')};
   def={'rs4886207'};
   dlgTitle='Input for SNP ID';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
	marker=char(answer);
	if (isempty(marker)),
	      return;
	elseif (length(marker)<4)
	      disp('Invalid Marker.')
	      return;
	elseif (~strcmp(marker([1 2]),'rs'))
	      disp('Invalid Marker.')
	      return;
	end
%marker='rs11571710';
%marker='rs6870660';
end



%figure;
if ~ishapmap3

    popset={'CEU','CHB','JPT','YRI'};
    nucset='ACGTN';
    for (k=1:4),
           [s] = snp_downloadhapmap(marker,char(popset{k}));
           filename=tempname;
           [fid,Msg] = fopen(filename,'wt');
           if fid == -1, error(Msg); end
           fprintf(fid,'%s',s);
           fclose(fid);
           [geno,genoinfo] = snp_readhapmap(filename);

           if isempty(geno), error('The marker is not in HapMap database, or genotype data is not avaiable!'); end
           %if (k==1|k==4),
           %[geno,genoinfo] = extrio(geno,genoinfo);
           %end
        [numHap,sizHap,seqHap] = counthaplotype(geno);
        sizHap=sizHap([1:min(3,length(sizHap))]);
        seqHap=seqHap([1:min(3,size(seqHap,1))],:);
            nuc=unique(seqHap);

        %fprintf('Population: %s  %2.2f %2.2f %2.2f\n', char(popset{k}),...
        %         sizHap(1)/sum(sizHap), sizHap(2)/sum(sizHap), sizHap(3)/sum(sizHap));


        genochar={};
        for (x=1:length(sizHap)),
              genochar{x}=char(nucset(seqHap(x,:)));
        end

        allechar={};
        allenum=zeros(1,length(nuc));
        for (y=1:length(nuc)),
              allenum(y)=sum(sum(geno==nuc(y)));
              allechar{y}=char(nucset(nuc(y)));
        end


        %subplot(2,4,k), pie([sizHap],{char(nucset(seqHap(1,:))),...
        %        char(nucset(seqHap(2,:))),char(nucset(seqHap(3,:)))});
        subplot(2,4,k), pie([sizHap], genochar)
        title(char(popset{k}))

        subplot(2,4,4+k),
        %pie([sum(sum(geno==nuc(1))),sum(sum(geno==nuc(2)))],allechar);
        pie(allenum, allechar);
            title(char(popset{k}))

    end

else
    %popset={'ASW','CEU','CHB','CHD','GIH','JPT','LWK','MEX','MKK','TSI','YRI'};
    popset={'LWK','MKK','YRI','ASW','CEU','TSI','GIH','CHB','CHD','JPT','MEX'};
    nucset='ACGTN';
    for (k=1:11),
           [s] = snp_downloadhapmap3(marker,popset{k});
           filename=tempname;
           [fid,Msg] = fopen(filename,'wt');
           if fid == -1, error(Msg); end
           fprintf(fid,'%s',s);
           fclose(fid);
           [geno,genoinfo] = snp_readhapmap(filename);

           if ~isempty(geno),

               %error('The marker is not in HapMap database, or genotype data is not avaiable!');

        [numHap,sizHap,seqHap] = counthaplotype(geno);
        sizHap=sizHap([1:min(3,length(sizHap))]);
        seqHap=seqHap([1:min(3,size(seqHap,1))],:);
            nuc=unique(seqHap);

        %fprintf('Population: %s  %2.2f %2.2f %2.2f\n', char(popset{k}),...
        %         sizHap(1)/sum(sizHap), sizHap(2)/sum(sizHap), sizHap(3)/sum(sizHap));


        genochar={};
        for (x=1:length(sizHap)),
              genochar{x}=char(nucset(seqHap(x,:)));
        end

        allechar={};
        allenum=zeros(1,length(nuc));
        for (y=1:length(nuc)),
              allenum(y)=sum(sum(geno==nuc(y)));
              allechar{y}=char(nucset(nuc(y)));
        end


        %subplot(2,4,k), pie([sizHap],{char(nucset(seqHap(1,:))),...
        %        char(nucset(seqHap(2,:))),char(nucset(seqHap(3,:)))});
        subplot(2,11,k), pie([sizHap], genochar)
        title(char(popset{k}))

        subplot(2,11,11+k),
        %pie([sum(sum(geno==nuc(1))),sum(sum(geno==nuc(2)))],allechar);
        pie(allenum, allechar);
          title(char(popset{k}))
           else
        subplot(2,11,k),
        title(char(popset{k}))
        subplot(2,11,11+k),
         title(char(popset{k}))
        end
    end


end

suptitle(strcat('SNP:',marker))

