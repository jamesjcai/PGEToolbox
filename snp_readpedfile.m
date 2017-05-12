function [genodata,markinfo] = snp_readpedfile(filename)

[genodata,markinfo]=snp_readlinkage(filename,'Delimiter','\t',...
            'MissingGenotype','0','UseACGT',true,'Noise',true);