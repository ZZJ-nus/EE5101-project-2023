c1T=C(1,:);
c2T=C(2,:);
%row1
num1=C(1,:)*B;%!=0
num2=C(1,:)*A*B;
index1=1;
%row2
num3=C(2,:)*B;%!=0
num4=C(2,:)*A*B;
index2=1;

B_star=[num1;num3];

C_star=[c1T*A;c2T*A];

C_star_star=[C(1,:)*(A+1*eye(4));C(2,:)*(A+1*eye(4))];
F=inv(B_star);
K=F*C_star_star;
[temp1,temp2]=eig(A-B*K);







