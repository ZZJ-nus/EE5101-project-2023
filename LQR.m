Q=eye(4);
R=[2 0;
    0 3];
T1=[A -1*B*inv(R)*transpose(B);
    -1*Q -1*transpose(A)];
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
K_LQR=inv(R)*transpose(B)*P;
%convert K into real number
K_LQR=real(K_LQR);




