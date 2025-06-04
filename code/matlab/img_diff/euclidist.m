function D=euclidist(data1,data2)
% EUCLIDIST: Euclidean distance between two points in 1-3 dimensions
%
%   Colour Engineering Toolbox
%   author:	© Phil Green
%   version:	1.1
%   date:	17-01-2001
%   book:	http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:     	http://www.digitalcolour.org

if size(data1)~=size(data2)
    error('Input coordinates have different dimensions');
end

c1=size(data1,2);
if c1==1
   D=abs(data1-data2);
elseif c1==2
   D=((data1(:,1)-data2(:,1)).^2+(data1(:,2)-data2(:,2)).^2).^0.5;
elseif c1==3
   D=((data1(:,1)-data2(:,1)).^2+(data1(:,2)-data2(:,2)).^2+(data1(:,3)-data2(:,3)).^2).^0.5;
else error('euclidist cannot calculate the distance with more than three coordinates')
end

