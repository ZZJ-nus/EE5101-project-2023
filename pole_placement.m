%multi input->one input
q=[1;2];
B1=B*q;
Wc=[B1 A*B1 A^2*B1 A^3*B1];
rank(Wc);
% Ackermann¡¯s formula
KT=[0 0 0 1]*inv(Wc)*poly(A,-0.8,-1);
K=q*KT;


%characteristic polynomial
function res=poly(A,pole1,pole2)
    res=(A^2 + 0.53*A + 0.19*eye(4))*(A - pole1*eye(4))*(A - pole2*eye(4));
end
    