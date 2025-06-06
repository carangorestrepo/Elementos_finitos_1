


clc
clear
syms xi eta x y

P =@(x,y)[[ 1, -x, -y, -x^2/2, - y^3/12 - (x^2*y)/2, -y^2/2, - x^3/12 - (y^2*x)/2, -(x*y)/2, x/2, (x*y)/2, y/2, (x*y)/2];
          [ 0,  1,  0,      x,                  x*y,      0,        x^2/4 + y^2/2,      y/2, 1/2,     y/2,   0,    -y/2];
          [ 0,  0,  1,      0,        x^2/2 + y^2/4,      y,                  x*y,      x/2,   0,    -x/2, 1/2,     x/2]];
      
Ci =@(xi,yi)[[1, -xi, -yi, -xi^2/2, -(xi^2*yi/2 + yi^3/12), -yi^2/2,-(xi*yi^2/2 + xi^3/12), -xi*yi/2, xi/2, xi*yi/2, yi/2, xi*yi/2];
             [0,   1,   0,      xi,                  xi*yi,       0,     (xi^2/4 + yi^2/2),     yi/2,  1/2,    yi/2,    0,   -yi/2];
             [0,   0,   1,       0,      (yi^2/4 + xi^2/2),      yi,                 xi*yi,     xi/2,    0,   -xi/2,  1/2,    xi/2]];
         
Q =@(eta,eta) [0, 0, 0, 1, eta, 0, xi/2, 0,  0,     0, 0,    0;
               0, 0, 0, 0,   0, 0,    0, 1, xi, eta/2, 0,    0;
               0, 0, 0, 0,   0, 0,    0, 0,  0,     0, 1, 2*xi;
               0, 0, 0, 0,   0, 0,    0, 0,  1,   eta, 0,    0;
               0, 0, 0, 0,   0, 0,    0, 0,  0,     0, 1,   xi];