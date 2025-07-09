function XYZ=lab2xyz_custom(Lab,XYZn)

% Lab2XYZ: Inverts the calculation of CIELAB L*a*b* 
% to find tristimulus XYZ from L*a*b* using the equation in
% Colorimetry and Colour Difference, in Green and MacDonald, (2002)
% Colour Engineering (Wiley).
%
% Example: XYZ=Lab2XYZ_k(Lab,XYZn)
% where XYZn is the reference white.
%
% If no argument is supplied for the reference white,
% D50 perfect diffuser is assumed. 
%
%   Colour Engineering Toolbox
%   author:    © Phil Green
%   version:   1.1
%   date:  	   16-02-2002
%   book:      http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:       http://www.digitalcolour.org

if nargin>1
   rwhite=XYZn;
else rwhite=[96.4212 100.0000 82.5188];
end

% set knee point of function
k=0.008856^(1/3);

L=Lab(:,1);a=Lab(:,2);b=Lab(:,3);
Xn=rwhite(1);Yn=rwhite(2);Zn=rwhite(3);

fy=(L+16)/116;
fx=a/500+fy;
fz=fy-b/200;

X=fx.^3*Xn;
Y=fy.^3*Yn;
Z=fz.^3*Zn;

% Now calculate T where T/Tn is below the knee point
p=fx<k;
q=fy<k;
r=fz<k;

X(p)=((fx(p)-16/116)/7.787)*Xn;
Y(q)=((fy(q)-16/116)/7.787)*Yn;
Z(r)=((fz(r)-16/116)/7.787)*Zn;

XYZ=[X,Y,Z];
