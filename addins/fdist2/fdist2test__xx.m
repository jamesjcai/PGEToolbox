load fdist2test geno12 geno1 geno2

%het2=snp_heterozygosity(geno2);
%het1=snp_heterozygosity(geno1);
%hom0=1-mean([het2,het1]);
%hom1=1-snp_heterozygosity(geno12);


[a,het2]=snp_heterozygosity(geno2);
[a,het1]=snp_heterozygosity(geno1);
[a,het12]=snp_heterozygosity(geno12);

h0=mean([het1;het2]);
h1=het12;


X=1-h0./h1;
Y=snp_fst(geno1,geno2);

figure;
plot(h1,Y,'o')


%{
loc_1	0.05	0.0	0.524874
loc_2	0.29	-0.043557	0.218571
loc_3	0.26	-0.012146	0.609548
loc_4	0.22	-0.04067	0.276462
loc_5	0.425	-0.027864	0.419009
loc_6	0.05	0.0	0.524874
loc_7	0.22	-0.04067	0.276462
loc_8	0.05	0.0	0.524874
loc_9	0.455556	-0.020841	0.484999
loc_10	0.425	-0.027864	0.419009
loc_11	0.425	-0.027864	0.419009
loc_12	0.05	0.0	0.524874
loc_13	0.18	-0.052632	0.068076
loc_14	0.275	0.100478	0.932008
loc_15	0.275	0.100478	0.932008
loc_16	0.18	-0.052632	0.068076
loc_17	0.185	0.004267	0.626761
loc_18	0.41	0.005135	0.679708
loc_19	0.41	0.005135	0.679708
loc_20	0.41	0.005135	0.679708
loc_21	0.275	0.100478	0.932008
loc_22	0.46	-0.029748	0.399399
loc_23	0.144444	-0.040486	0.279065
loc_24	0.18	-0.052632	0.068076
loc_25	0.18	-0.052632	0.068076
loc_26	0.18	-0.052632	0.068076
loc_27	0.46	-0.029748	0.399399
loc_28	0.41	0.005135	0.679708
loc_29	0.111111	0.058824	0.829778
loc_30	0.4375	-0.028571	0.41174
loc_31	0.1	0.052632	0.812631
loc_32	0.46	-0.029748	0.399399
loc_33	0.05	0.0	0.524874
loc_34	0.05	0.0	0.524874
loc_35	0.425	-0.05934	0.063596
loc_36	0.41	0.005135	0.679708
loc_37	0.275	0.100478	0.932008
loc_38	0.47	-0.047032	0.195
loc_39	0.05	0.0	0.524874
loc_40	0.41	0.005135	0.679708
loc_41	0.46	-0.029748	0.399399
loc_42	0.302469	0.114046	0.94424

%}