function [w]=Q3dBasis(x,y,index)

q1=.5-.5*sqrt((3+2*sqrt(6/5))/7);
q2=.5-.5*sqrt((3-2*sqrt(6/5))/7);
q3=.5+.5*sqrt((3-2*sqrt(6/5))/7);
q4=.5+.5*sqrt((3+2*sqrt(6/5))/7);

w=[ (x-q2)*(x-q3)*(x-q4)*(y-q2)*(y-q3)*(y-q4)/((q1-q2)*(q1-q3)*(q1-q4)*(q1-q2)*(q1-q3)*(q1-q4))
    (x-q1)*(x-q3)*(x-q4)*(y-q2)*(y-q3)*(y-q4)/((q2-q1)*(q2-q3)*(q2-q4)*(q1-q2)*(q1-q3)*(q1-q4))
    (x-q1)*(x-q2)*(x-q4)*(y-q2)*(y-q3)*(y-q4)/((q3-q1)*(q3-q2)*(q3-q4)*(q1-q2)*(q1-q3)*(q1-q4))
    (x-q1)*(x-q2)*(x-q3)*(y-q2)*(y-q3)*(y-q4)/((q4-q1)*(q4-q2)*(q4-q3)*(q1-q2)*(q1-q3)*(q1-q4))
    (x-q2)*(x-q3)*(x-q4)*(y-q1)*(y-q3)*(y-q4)/((q1-q2)*(q1-q3)*(q1-q4)*(q2-q1)*(q2-q3)*(q2-q4))
    (x-q1)*(x-q3)*(x-q4)*(y-q1)*(y-q3)*(y-q4)/((q2-q1)*(q2-q3)*(q2-q4)*(q2-q1)*(q2-q3)*(q2-q4))
    (x-q1)*(x-q2)*(x-q4)*(y-q1)*(y-q3)*(y-q4)/((q3-q1)*(q3-q2)*(q3-q4)*(q2-q1)*(q2-q3)*(q2-q4))
    (x-q1)*(x-q2)*(x-q3)*(y-q1)*(y-q3)*(y-q4)/((q4-q1)*(q4-q2)*(q4-q3)*(q2-q1)*(q2-q3)*(q2-q4))
    (x-q2)*(x-q3)*(x-q4)*(y-q1)*(y-q2)*(y-q4)/((q1-q2)*(q1-q3)*(q1-q4)*(q3-q1)*(q3-q2)*(q3-q4))
    (x-q1)*(x-q3)*(x-q4)*(y-q1)*(y-q2)*(y-q4)/((q2-q1)*(q2-q3)*(q2-q4)*(q3-q1)*(q3-q2)*(q3-q4))
    (x-q1)*(x-q2)*(x-q4)*(y-q1)*(y-q2)*(y-q4)/((q3-q1)*(q3-q2)*(q3-q4)*(q3-q1)*(q3-q2)*(q3-q4))
    (x-q1)*(x-q2)*(x-q3)*(y-q1)*(y-q2)*(y-q4)/((q4-q1)*(q4-q2)*(q4-q3)*(q3-q1)*(q3-q2)*(q3-q4))
    (x-q2)*(x-q3)*(x-q4)*(y-q1)*(y-q2)*(y-q3)/((q1-q2)*(q1-q3)*(q1-q4)*(q4-q1)*(q4-q2)*(q4-q3))
    (x-q1)*(x-q3)*(x-q4)*(y-q1)*(y-q2)*(y-q3)/((q2-q1)*(q2-q3)*(q2-q4)*(q4-q1)*(q4-q2)*(q4-q3))
    (x-q1)*(x-q2)*(x-q4)*(y-q1)*(y-q2)*(y-q3)/((q3-q1)*(q3-q2)*(q3-q4)*(q4-q1)*(q4-q2)*(q4-q3))
    (x-q1)*(x-q2)*(x-q3)*(y-q1)*(y-q2)*(y-q3)/((q4-q1)*(q4-q2)*(q4-q3)*(q4-q1)*(q4-q2)*(q4-q3))];
    
    