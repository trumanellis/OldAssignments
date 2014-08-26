function w=Q1dBasis(x,y,index)

w=[ 3*(x-.5-.5/sqrt(3))*(y-.5-.5/sqrt(3));
    -3*(x-.5+.5/sqrt(3))*(y-.5-.5/sqrt(3));
    3*(x-.5+.5/sqrt(3))*(y-.5+.5/sqrt(3));
    -3*(x-.5-.5/sqrt(3))*(y-.5+.5/sqrt(3))];
w(abs(w) < 1e-10)=0;