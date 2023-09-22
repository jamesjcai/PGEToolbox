function [ancalle, devalle] = snp_ancalleweb(rsid)
%SNP_ANCALLEWEB - returns major allele and frequency for ancestral GENO

ancalle = [];
devalle = [];
[chrid, pos] = snp_locator(rsid);
if chrid > 0
    ancalle = snp_ancallechimp(chrid, pos);
    devalle = encodeseq(chrnuc(chrid, pos));
end
