X_SP=[0;0.5;-0.4;0.3];
Y_SP=C*X_SP;

%argument system 
%x_argument 6*1
A_aug=[A zeros(4,2);
    -1*C zeros(2,2)];%6*6
B_aug=[B;zeros(2,2)];%6*2
C_aug=[C zeros(2,2)];%2*6
QT=[A B;C zeros(2,2)];

%%%LQR  QX RU
Q=eye(6);
R=[0.02 0;
    0 0.03];
T1=[A_aug -1*B_aug*inv(R)*transpose(B_aug);
    -1*Q -1*transpose(A_aug)];
[V1,D1]=eig(T1);% eig vector+eig val
val=diag(D1);%eig val
[row,col]=size(V1);
% initialize 2 empty vector
temp1=[];
temp2=[];
for i=1:col %each col is a eig vector
    if(val(i,:)<0) % eig val must be stable
        temp=V1(:,i); % eig vector
        temp1=horzcat(temp1,temp(row/2+1:row,:));
        temp2=horzcat(temp2,temp(1:row/2,:));
    end
end
P=temp1*inv(temp2);
K_LQR=inv(R)*transpose(B_aug)*P;
%convert K into real number
K_LQR=real(K_LQR);
K1=K_LQR(:,1:4);
K2=K_LQR(:,5:6);

