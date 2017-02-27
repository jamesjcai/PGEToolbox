function [d]=braycurtisdist(v1,v2)
%Bray Curtis distance
%http://people.revoledu.com/kardi/tutorial/Similarity/index.html

a=sum(sum(abs(v1-v2)));
b=sum(sum(v1+v2));
d=a./b;

%Bray Curtis distance sometimes is also called Sorensen distance is a
%normalization method that common used in botany, ecology and environmental
%sicence field. It views the space as grid similar to the city block
%distance. The Bray curtis distance has a nice property that if all
%coordinates are postive, its value is between zero and one. Zero bray
%curtis represent exact similar coordinate. If both objects are in the zero
%coordinates, the Bray curtis distance is undefined. The normalization is
%done using absolute difference divided by the summation.








