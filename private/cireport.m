function [a,b]=cireport(x)
a=prctile(x,2.5);
b=prctile(x,97.5);
disp(['Mean:                                ',num2str(mean(x))]);
disp(['Standard Deviation:                  ',num2str(std(x))]);
disp(['Median:                              ',num2str(median(x))]);
disp(['2.5th Percentile:                     ',num2str(a)]);
disp(['97.5th Percentile:                     ',num2str(b)]);