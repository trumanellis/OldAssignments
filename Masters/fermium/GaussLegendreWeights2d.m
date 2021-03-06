function [pts2d,wgts2d] = GaussLegendreWeights2d(order)
% This actually uses Gauss-Legendre for order 3 and higher

switch order
    case 1
        pts2d=[.5 .5];
        wgts2d=1;
    case 2
        pts2d=[.5-.5/sqrt(3) .5-.5/sqrt(3)
               .5+.5/sqrt(3) .5-.5/sqrt(3)
               .5+.5/sqrt(3) .5+.5/sqrt(3)
               .5-.5/sqrt(3) .5+.5/sqrt(3)];
%         pts2d=[ 0.21132487 0.21132487
%                 0.21132487 0.78867513
%                 0.78867513 0.21132487
%                 0.78867513 0.78867513];

        wgts2d=[.25; .25; .25; .25];
    case 3
        pts2d=[.5-.5*sqrt(3/5) .5-.5*sqrt(3/5)
               .5+.5*sqrt(3/5) .5-.5*sqrt(3/5)
               .5+.5*sqrt(3/5) .5+.5*sqrt(3/5)
               .5-.5*sqrt(3/5) .5+.5*sqrt(3/5)
               .5              .5-.5*sqrt(3/5)
               .5+.5*sqrt(3/5) .5
               .5              .5+.5*sqrt(3/5)
               .5-.5*sqrt(3/5) .5
               .5              .5];
%         pts2d=[ 0.11270167 0.11270167
%                 0.11270167 0.5
%                 0.11270167 0.88729833
%                 0.5 0.11270167
%                 0.5 0.5
%                 0.5 0.88729833
%                 0.88729833 0.11270167
%                 0.88729833 0.5
%                 0.88729833 0.88729833];
            wgts2d=[ 25/324  
                     25/324
                     25/324
                     25/324
                     10/81  
                     10/81  
                     10/81  
                     10/81
                     16/81];
%              wgts2d=[ 0.07716049
%                      0.12345679  
%                      0.07716049  
%                      0.12345679  
%                      0.19753086  
%                      0.12345679
%                      0.07716049  
%                      0.12345679  
%                      0.07716049];
    case 4
        pts2d=[.5-.5*sqrt((3+2*sqrt(6/5))/7) .5-.5*sqrt((3+2*sqrt(6/5))/7)
               .5-.5*sqrt((3-2*sqrt(6/5))/7) .5-.5*sqrt((3+2*sqrt(6/5))/7)
               .5+.5*sqrt((3-2*sqrt(6/5))/7) .5-.5*sqrt((3+2*sqrt(6/5))/7)
               .5+.5*sqrt((3+2*sqrt(6/5))/7) .5-.5*sqrt((3+2*sqrt(6/5))/7)
               .5-.5*sqrt((3+2*sqrt(6/5))/7) .5-.5*sqrt((3-2*sqrt(6/5))/7)
               .5-.5*sqrt((3-2*sqrt(6/5))/7) .5-.5*sqrt((3-2*sqrt(6/5))/7)
               .5+.5*sqrt((3-2*sqrt(6/5))/7) .5-.5*sqrt((3-2*sqrt(6/5))/7)
               .5+.5*sqrt((3+2*sqrt(6/5))/7) .5-.5*sqrt((3-2*sqrt(6/5))/7)
               .5-.5*sqrt((3+2*sqrt(6/5))/7) .5+.5*sqrt((3-2*sqrt(6/5))/7)
               .5-.5*sqrt((3-2*sqrt(6/5))/7) .5+.5*sqrt((3-2*sqrt(6/5))/7)
               .5+.5*sqrt((3-2*sqrt(6/5))/7) .5+.5*sqrt((3-2*sqrt(6/5))/7)
               .5+.5*sqrt((3+2*sqrt(6/5))/7) .5+.5*sqrt((3-2*sqrt(6/5))/7)
               .5-.5*sqrt((3+2*sqrt(6/5))/7) .5+.5*sqrt((3+2*sqrt(6/5))/7)
               .5-.5*sqrt((3-2*sqrt(6/5))/7) .5+.5*sqrt((3+2*sqrt(6/5))/7)
               .5+.5*sqrt((3-2*sqrt(6/5))/7) .5+.5*sqrt((3+2*sqrt(6/5))/7)
               .5+.5*sqrt((3+2*sqrt(6/5))/7) .5+.5*sqrt((3+2*sqrt(6/5))/7)];
           wgts2d=[ (354-36*sqrt(30))/5184
                    294/5184
                    294/5184
                    (354-36*sqrt(30))/5184
                    294/5184
                    (354+36*sqrt(30))/5184
                    (354+36*sqrt(30))/5184
                    294/5184
                    294/5184
                    (354+36*sqrt(30))/5184
                    (354+36*sqrt(30))/5184
                    294/5184
                    (354-36*sqrt(30))/5184
                    294/5184
                    294/5184
                    (354-36*sqrt(30))/5184];
%         pts2d=[ 0.06943184 0.06943184  
%                 0.06943184 0.33000948  
%                 0.06943184 0.66999052 
%                 0.06943184 0.93056816
%                 0.33000948 0.06943184  
%                 0.33000948 0.33000948  
%                 0.33000948 0.66999052 
%                 0.33000948 0.93056816
%                 0.66999052 0.06943184  
%                 0.66999052 0.33000948  
%                 0.66999052 0.66999052 
%                 0.66999052 0.93056816
%                 0.93056816 0.06943184  
%                 0.93056816 0.33000948  
%                 0.93056816 0.66999052 
%                 0.93056816 0.93056816];
%             wgts2d=[ 0.10632333  
%                      0.05671296  
%                      0.05671296  
%                      0.10632333  
%                      0.05671296  
%                      0.03025075
%                      0.03025075  
%                      0.05671296  
%                      0.05671296  
%                      0.03025075  
%                      0.03025075  
%                      0.05671296
%                      0.10632333  
%                      0.05671296  
%                      0.05671296  
%                      0.10632333];
    case 5
        pts2d=[ 0.04691008 0.04691008
                0.04691008 0.23076534 
                0.04691008 0.5 
                0.04691008 0.76923466
                0.04691008 0.95308992
                0.23076534 0.04691008
                0.23076534 0.23076534 
                0.23076534 0.5 
                0.23076534 0.76923466
                0.23076534 0.95308992
                0.5  0.04691008
                0.5  0.23076534 
                0.5  0.5 
                0.5  0.76923466
                0.5  0.95308992
                0.76923466 0.04691008
                0.76923466 0.23076534 
                0.76923466 0.5 
                0.76923466 0.76923466
                0.76923466 0.95308992
                0.95308992 0.04691008
                0.95308992 0.23076534 
                0.95308992 0.5 
                0.95308992 0.76923466
                0.95308992 0.95308992];
            wgts2d=[ 0.01403359  
                     0.02835
                     0.03369627
                     0.02835  
                     0.01403359
                     0.02835
                     0.05727135 
                     0.06807163
                     0.05727135 
                     0.02835
                     0.03369627  
                     0.06807163
                     0.08090864  
                     0.06807163
                     0.03369627
                     0.02835
                     0.05727135
                     0.06807163
                     0.05727135
                     0.02835 
                     0.01403359 
                     0.02835 
                     0.03369627 
                     0.02835
                     0.01403359];

%         pts2d=[
%                          0                         0
%                          0         0.172673164646011
%                          0                       0.5
%                          0         0.827326835353989
%                          0                         1
%          0.172673164646011                         0
%          0.172673164646011         0.172673164646011
%          0.172673164646011                       0.5
%          0.172673164646011         0.827326835353989
%          0.172673164646011                         1
%                        0.5                         0
%                        0.5         0.172673164646011
%                        0.5                       0.5
%                        0.5         0.827326835353989
%                        0.5                         1
%          0.827326835353989                         0
%          0.827326835353989         0.172673164646011
%          0.827326835353989                       0.5
%          0.827326835353989         0.827326835353989
%          0.827326835353989                         1
%                          1                         0
%                          1         0.172673164646011
%                          1                       0.5
%                          1         0.827326835353989
%                          1                         1];
%          wgts2d=[
%                     0.0025
%         0.0136111111111111
%         0.0177777777777778
%         0.0136111111111111
%                     0.0025
%         0.0136111111111111
%         0.0741049382716049
%         0.0967901234567901
%         0.0741049382716049
%         0.0136111111111111
%         0.0177777777777778
%         0.0967901234567901
%           0.12641975308642
%         0.0967901234567901
%         0.0177777777777778
%         0.0136111111111111
%         0.0741049382716049
%         0.0967901234567901
%         0.0741049382716049
%         0.0136111111111111
%                     0.0025
%         0.0136111111111111
%         0.0177777777777778
%         0.0136111111111111
%                     0.0025];
end