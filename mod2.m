%A code for analysis of laminate
clc;
clear all;
prompt1="Enter the no. of layers : ";
 n=input(prompt1);
%disp('Enter the value of E_1,E_2,miu12,G_12 for each layer');
coe=zeros((4*n),1);
%for i=1:4*n 
 %   coe(i)=input("constant : ");
%end
coe=[181*10^9,10.3*10^9,0.28,7.17*10^9,181*10^9,10.3*10^9,0.28,7.17*10^9,181*10^9,10.3*10^9,0.28,7.17*10^9];
strength=zeros((5*n),1);
%disp('Enter the value of sigma_1_Tu,sigma_1_Cu,sigma_2_Tu,sigma_2_Cu,Tau_12_u respectively for each layer');
%for i=1:5*n 
 %   strength(i)=input("constant : ");
%end
strength=[1500*10^6,1500*10^6,40*10^6,246*10^6,68*10^6,1500*10^6,1500*10^6,40*10^6,246*10^6,68*10^6,1500*10^6,1500*10^6,40*10^6,246*10^6,68*10^6];
%Reduced stiffness matrix Q
Q=zeros(4*n,1);Qm=zeros(3*n,3);
countq=1;
for i=1:n
    d=1-coe(3+4*(i-1))*(coe(3+4*(i-1))*coe(2+4*(i-1))/coe(1+4*(i-1)));
    Q(countq) = coe(1+4*(i-1))/d;
    Qm(1+3*(i-1),1)=Q(countq);
    countq=countq+1;
    Q(countq)=coe(2+4*(i-1))*coe(3+4*(i-1))/d;
    Qm(1+3*(i-1),2)=Q(countq);
    Qm(2+3*(i-1),1)=Q(countq);
    countq=countq+1;
    Q(countq)=coe(2+4*(i-1))/d;
    Qm(2+3*(i-1),2)=Q(countq);
    countq=countq+1;
    Q(countq)=coe(4+4*(i-1))
    Qm(3+3*(i-1),3)=Q(countq);
    countq=countq+1;
end
Qmod=Qm;
%constructing Qm for odd layers
if(rem(n,2)~=0)
    Qmod=zeros(3*(n+1),3);
    for i=1:n
        if(i==ceil(n/2))
            Qmod(1+3*(i-1):3+3*(i-1),:)=Qm(1+3*(i-1):3+3*(i-1),:);
            Qmod(1+3*(i):3+3*(i),:)=Qm(1+3*(i-1):3+3*(i-1),:);
            
        elseif(i>ceil(n/2))
            Qmod(1+3*(i):3+3*(i),:)=Qm(1+3*(i-1):3+3*(i-1),:);
        else
            Qmod(1+3*(i-1):3+3*(i-1),:)=Qm(1+3*(i-1):3+3*(i-1),:);
        end
    end
end
disp("Enter the orientation angle in degree and thickness of each layer");
theta=zeros(n,1);
th=zeros(n,1);
for i=1:n
    theta(i)=input("angle value : ");
    th(i)=input("thickness value : ");
end
%Q_bar
countqbar=1;
Q_bar=zeros(6*n,1);
for i=1:n
    s=sind(theta(i));
    c=cosd(theta(i));
Q_bar(countqbar)=Q(1+4*(i-1))*c^4 + 2*(Q(2+4*(i-1)) + 2*Q(4+4*(i-1)))*(c*s)^2 + Q(3+4*(i-1))*s^4; %Q11
countqbar=countqbar+1;
Q_bar(countqbar)=(Q(1+4*(i-1)) +Q(3+4*(i-1)) - 4*Q(4+4*(i-1)))*(c*s)^2 + Q(2+4*(i-1))*(s^4 + c^4); %Q12
countqbar=countqbar+1;
Q_bar(countqbar)=(Q(1+4*(i-1)) -Q(2+4*(i-1)) -2*Q(4+4*(i-1)))*(s*c^3) - (Q(3+4*(i-1)) -Q(2+4*(i-1)) -2*Q(4+4*(i-1)))*(c*s^3); %Q16
countqbar=countqbar+1;
Q_bar(countqbar)=Q(1+4*(i-1))*s^4 + 2*(Q(2+4*(i-1)) + 2*Q(4+4*(i-1)))*(c*s)^2 + Q(3+4*(i-1))*c^4; %Q22
countqbar=countqbar+1;
Q_bar(countqbar)=(Q(1+4*(i-1)) -Q(2+4*(i-1)) -2*Q(4+4*(i-1)))*(c*s^3) - (Q(3+4*(i-1)) -Q(2+4*(i-1)) -2*Q(4+4*(i-1)))*(s*c^3); %Q26
countqbar=countqbar+1;
Q_bar(countqbar)=(Q(1+4*(i-1)) + Q(3+4*(i-1)) - 2*Q(2+4*(i-1)) - 2*Q(4+4*(i-1)))*(c*s)^2 + Q(4+4*(i-1))*(c^4 + s^4) %Q66
countqbar=countqbar+1;
end
if(rem(n,2)==0)
    zt=zeros(n/2,1);
    zb=zeros(n/2,1);
    st=0;sb=0;
    for i=1:n/2
        zb(i)=th(n/2 + i);
        sb=sb+zb(i);
        zb(i)=sb;
        zt(i)= -th(n/2 -(i-1));
        st=st+zt(i);
        zt(i)=st;
    end
else
    zt=zeros(fix(n/2 + 1),1);
    zb=zeros(fix(n/2 + 1),1);
    st=0;sb=0,counter=1;
    for i=1:fix(n/2 +1)
        if (i==1)
            zb(i)=th(fix(n/2 +1))/2;
            sb=sb+zb(i);
            zt(i)=-th(fix(n/2 +1))/2;
            st=st+zt(i);
        else
            zb(i)=th(ceil(n/2) + counter);
            sb=sb+zb(i);
            zb(i)=sb;
            zt(i)= -th(ceil(n/2) -(counter));
            st=st+zt(i);
            zt(i)=st;
            counter=counter+1;
        end
    end
end
cn=1;
Q_b=zeros(3*n,3);
%constructing Q_b from vector Q_bar
for i=1:n
    for k=1:3
        for l=1:3
            if(l>=k)
            Q_b(k+(3*(i-1)),l)=Q_bar(cn);
            cn=cn+1;
            if(k~=l)
                Q_b(l+(3*(i-1)),k)=Q_b(k+(3*(i-1)),l);
            end
            end
        end
    end
end
Q_bmod=Q_b;
if(rem(n,2)~=0)
    Q_bmod=zeros(3*(n+1),3);
    for i=1:n
        if(i==ceil(n/2))
            Q_bmod(1+3*(i-1):3+3*(i-1),:)=Q_b(1+3*(i-1):3+3*(i-1),:);
            Q_bmod(1+3*(i):3+3*(i),:)=Q_b(1+3*(i-1):3+3*(i-1),:);
            
        elseif(i>ceil(n/2))
            Q_bmod(1+3*(i):3+3*(i),:)=Q_b(1+3*(i-1):3+3*(i-1),:);
        else
            Q_bmod(1+3*(i-1):3+3*(i-1),:)=Q_b(1+3*(i-1):3+3*(i-1),:);
        end
    end
end
A=zeros(3,3);
B=zeros(3,3);D=zeros(3,3);
zt=flip(zt);
z=[zt;zb];lenz=length(z);
if (rem(n,2)==0)
    z=[zt;0;zb];
end
for i=1:n
        A= A+ Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)-z(i)); 
        B= B+ Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)^2-z(i)^2);
        D= D+ Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)^3-z(i)^3);
end
B=B./2;
D=D./3;
R=[A B;B D];
NM=[100000 0 0 0 0 0];
sr=inv(R)*NM'; %mid strain and curvature vector
midst=sr(1:3);
K=sr(4:6);
if(rem(n,2)~=0)
    thet=zeros(length(z),1);
    for i=1:length(z)
        if(i==ceil(n/2))
            thet(i)=theta(ceil(n/2));
            thet(i+1)=theta(ceil(n/2));
            i=i+1;
        elseif(i>ceil(n/2))
            thet(i)=theta(i-1);
        else
            thet(i)=theta(i);
        end
    end
end
epxy=zeros(3*lenz,1);
T=zeros(3*lenz,3);   %transformation matrix
for i=1:lenz
    epxy(1+3*(i-1):3+3*(i-1)) = midst + z(i)*K;
    epxy(3*i)=epxy(3*i)/2;    %making gammaxy to gammaxy/2
    if (rem(n,2)==0)
    c=cosd(theta(i));s=sind(theta(i));
    else
        c=cosd(thet(i));s=sind(thet(i));
    end
    %defining elements of trans formation matrix
    T(1+3*(i-1),1)=c^2;
    T(1+3*(i-1),2)=s^2;
    T(1+3*(i-1),3)=2*s*c;
    T(2+3*(i-1),1)=s^2;
    T(2+3*(i-1),2)=c^2;
    T(2+3*(i-1),3)=-2*s*c;
    T(3+3*(i-1),1)=-s*c;
    T(3+3*(i-1),2)=s*c;
    T(3+3*(i-1),3)=c^2-s^2;
end
ep12=zeros(3*lenz,1);
for i=1:lenz
    ep12(1+3*(i-1):3+3*(i-1))=T(1+3*(i-1):3+3*(i-1),:)*epxy(1+3*(i-1):3+3*(i-1));
    ep12(3*i)=ep12(3*i)*2;    %gamma12/2 to gamma12
end
sigma12=zeros(3*lenz,1);
for i=1:lenz
    sigma12(1+3*(i-1):3+3*(i-1))=Qmod(1+3*(i-1):3+3*(i-1),:)*ep12(1+3*(i-1):3+3*(i-1));
end
sratio=zeros(3*lenz,1);    %Finding strength ratios to apply maximum stress theory
for i=1:n
    if (sigma12(i+3*(i-1))>=0)
        sratio(i+3*(i-1))=sigma12(i+3*(i-1))/strength(i+3*(i-1));
    else
        sratio(i+3*(i-1))=-sigma12(i+3*(i-1))/strength(i+1+3*(i-1));
    end
    if (sigma12(i+1+3*(i-1))>=0)
        sratio(i+1+3*(i-1))=sigma12(i+1+3*(i-1))/strength(i+2+3*(i-1));
    else
        sratio(i+1+3*(i-1))=-sigma12(i+1+3*(i-1))/strength(i+3+3*(i-1));
    end
    sratio(i+2+3*(i-1))=sigma12(i+2+3*(i-1))/strength(i+4+3*(i-1));
end
maxsr=sratio(1);maxsr_k=1;
for i=1:3*n
    if(sratio(i)>maxsr)
        maxsr=sratio(i);maxsr_k=i;
    end
end
FPF=NM./maxsr
Failed_layer=ceil(maxsr_k/3)
failed_theta=theta(Failed_layer)
failed_layers=(theta==failed_theta)
B=B.*2;
D=D.*3;
for i=1:n
    if (failed_layers(i)==1)
        A= A - Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)-z(i));
        B= B - Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)^2-z(i)^2);
        D= D - Q_b(1+3*(i-1):3+3*(i-1),:).*(z(i+1)^3-z(i)^3);
        Q_b(1+3*(i-1):3+3*(i-1),:)=zeros(3,3);
    end
end
B=B./2;
D=D./3;
R=[A B;B D];
sr=inv(R)*FPF'; %mid strain and curvature vector
midst=sr(1:3);
K=sr(4:6);
for i=1:lenz
    epxy(1+3*(i-1):3+3*(i-1)) = midst + z(i)*K;
    epxy(3*i)=epxy(3*i)/2;    %making gammaxy to gammaxy/2
end

