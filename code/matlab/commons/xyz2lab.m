function Lab=xyz2lab(XYZ,XYZn)

% XYZ2LAB: calculates L,a,b from XYZ according to CIE 15:2004
%
% Example: LAB=XYZ2Lab(XYZ,XYZn)
% where XYZn is the reference white.
%
% If no argument is supplied for the reference white, the default 
% D50 is used.
% Note that by calling the d.m function, the reference white can be
% specified as follows:
%   LAB=XYZ2Lab(XYZ,d(50))  % CIE illuminant D50
%   LAB=XYZ2Lab(XYZ,d('A')) % CIE illuminant A
% > help d for details of illuminants supported.
%
%   Colour Engineering Toolbox
%   author:    © Phil Green
%   version:   1.2
%   date:  	   16-02-2002 (CIE 15.2 equations)
%   revised    19-09-2007 (CIE 15:2004 equations)
%   book:      http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:       http://www.digitalcolour.org

if nargin<2
   rwhite=d(50);
else
   rwhite=XYZn;
end

Xn=rwhite(1);Yn=rwhite(2);Zn=rwhite(3); 
X=XYZ(:,1);Y=XYZ(:,2);Z=XYZ(:,3);

% precompute constant for knee function
const=(24/116)^3;

% normalize to reference white
Yrel=Y/Yn;
Xrel=X/Xn;
Zrel=Z/Zn;

% apply cube root
fY=Yrel.^(1/3);
fX=Xrel.^(1/3);
fZ=Zrel.^(1/3);

% find T/Tn <= 24/116
p=find(Yrel<=const);
q=find(Xrel<=const);
r=find(Zrel<=const);

% calculate f(T) for T/Tn <= (24/116)^3
fY(p)=(841/108)*Yrel(p)+16/116;
fX(q)=(841/108)*Xrel(q)+16/116;
fZ(r)=(841/108)*Zrel(r)+16/116;

% calculate L,a,b
L=116*fY-16;
a=500*(fX-fY);
b=200*(fY-fZ);

Lab=[L,a,b];
