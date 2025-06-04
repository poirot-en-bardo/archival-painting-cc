function DE94=cie94(LABREF,LAB,K)

%   CIE94: Computes CIE94 colour difference between a reference and a sample
%   cie94.m accepts:
%   -  a single reference and single sample value
%   -  multiple values for both reference and sample
%   -  single reference and multiple samples
%
%   The weighting parameters kL,kC,kH can be supplied as optional arguments;
%   otherwise they default to 1
%
%   Colour Engineering Toolbox
%   author:    © Phil Green
%   version:   1.2
%   date:  	   10-07-2002
%   book:      http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:       http://www.digitalcolour.org

if nargin>2
   if length(K)>2
      kL=K(1);kC=K(2);kH=K(3);
   end
else
   kL=1;kC=1;kH=1;
end

Lref=LABREF(:,1);aref=LABREF(:,2);bref=LABREF(:,3);
Cref=(aref.^2+bref.^2).^0.5;

L=LAB(:,1);a=LAB(:,2);b=LAB(:,3);
C=(a.^2+b.^2).^0.5;

Sc=1+0.045*(Cref+C)/2;
Sh=1+0.015*(Cref+C)/2;

DC=abs(Cref-C);
DL=abs(Lref-L);
DE=(DL.^2+(aref-a).^2+(bref-b).^2).^0.5;
DH=real((DE.^2-DL.^2-DC.^2).^0.5);

DE94=((DL/kL).^2+(DC./(kC*Sc)).^2+(DH./(kH*Sh)).^2).^0.5;