function [AREA,XCEN,YCEN]=propiedadesgeo(xy)
AREA = 0;
XCEN = 0;
YCEN = 0;
%IXX = 0;
%IYY = 0;
%IXY = 0;
N = size(xy(:,1));
for I = 1: N - 1
    X1 = xy(I, 1);
    Y1 = xy(I, 2);
    X2 = xy(I + 1, 1);
    Y2 = xy(I + 1, 2);
    AREA = AREA + (Y2 - Y1) * (X2 + X1) / 2;
    XCEN = XCEN + (Y2 - Y1) / 8 * ((X2 + X1) ^ 2 + (X2 - X1) ^ 2 / 3);
    YCEN = YCEN + (X2 - X1) / 8 * ((Y2 + Y1) ^ 2 + (Y2 - Y1) ^ 2 / 3);
    %IXX = IXX + (X2 - X1) * (Y2 + Y1) / 24 * ((Y2 + Y1) ^ 2 + (Y2 - Y1) ^ 2);
    %IYY = IYY + (Y2 - Y1) * (X2 + X1) / 24 * ((X2 + X1) ^ 2 + (X2 - X1) ^ 2);
    %DIFER = X2 - X1;
    %if DIFER == 0 
    %    DIFER = 0.000001;        
    %end
%    IXY = IXY + (1 / 8 * (Y2 - Y1) ^ 2 * (X2 + X1) * (X2 ^ 2 + X1 ^ 2) + 1 / 3 * (Y2 - Y1) * (X2 * Y1 - X1 * Y2) * (X2 ^ 2 + X2 * X1 + X1 ^ 2) + 1 / 4 * (X2 * Y1 - X1 * Y2) ^ 2 * (X2 + X1)) / DIFER;
end
AREA = -AREA;
XCEN = -XCEN / AREA;
YCEN = YCEN / AREA;
%IYY = -IYY;
%IXXC = IXX - AREA * YCEN ^ 2;
%IYYC = IYY - AREA * XCEN ^ 2;
%IXYC = IXY - AREA * XCEN * YCEN;
%RESTA = IXXC - IYYC;
%if RESTA == 0 
%    RESTA = 0.000001;
%end 
%TETA = 0.5 * atan(-2 * IXYC / RESTA);
%Pi = 4 * atan(1);
%TETA = TETA * 180 / Pi;
%R = (IXXC / AREA) ^ (1 / 2);