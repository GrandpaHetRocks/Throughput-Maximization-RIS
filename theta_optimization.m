clc
clear all
% rng(65)
k=4;
n=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=randn(1,1);
% b=randn(1,1);
% hdk=complex(a,b);
% n=10;
% a=randn(1,n);
% b=randn(1,n);
% theta=complex(a,b)';
% a=randn(n,1);
% b=randn(n,1);
% Hrk=complex(a,b);
% hk=hdk+theta*Hrk
%%%%%%%%%%%%%%for one client%%%%%%%%%%
P=ones(1,k)*100;
a=randn(k,1);
b=randn(k,1);
hdk_array=complex(a,b);
a=randn(1,n);
b=randn(1,n);
theta=complex(a,b)';
theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));
Hrk_array=[];
for i=1:k
a=randn(n,1);
b=randn(n,1);
Hrk=complex(a,b);
Hrk_array=[Hrk_array Hrk];
end

alphak=0;
sumrate=[];
for i=1:50
    sum3=0;
    for j=2:k
        sum3=sum3+abs(hdk_array(j)'+theta'*Hrk_array(:,j))^2*P(j);
    end
    gamma=abs(hdk_array(1)'+theta'*Hrk_array(:,1))^2*P(1)/sum3;
    sumrate=[sumrate log2(1+gamma)];
    sum=0;
    
    for j=1:k    
        hk_=hdk_array(j)'+theta'*Hrk_array(:,j)*P(j);
        hk_=abs(hk_)^2;
        sum=sum+hk_;
    end
    betak=(sqrt(P(1)*(1+alphak)) * (hdk_array(1)' + theta'*Hrk_array(:,1))*P(1))/sum;
    hk=hdk_array(1)+Hrk_array(:,1)'*theta;
    kaik=1/sqrt(P(1))*real(conj(betak)*hk'*P(1));
    alphak=(kaik^2+kaik*sqrt(kaik^2+4))/2;
    sum1=0;
    for j=1:k
        a1=Hrk_array(:,j)*P(j);
        sum1=sum1+a1'*a1;
    end
    
    U=abs(betak)^2*sum1;
    
    sum2=0;
    for j=1:k
        b1=hdk_array(j)*P(j);
        a1=Hrk_array(:,j)*P(j);
        sum2=sum2+ conj(b1)*a;
    end
    
    v=sqrt(P(1)*(1+alphak))*conj(betak)*Hrk_array(:,1)*P(1) - abs(betak')^2*sum2;
    
    phi=angle(theta);
    temp1=U*exp(j*phi)-v;
    grad=2*real(-j*exp(-j*temp1));
    phi=phi-grad/10000000;
    theta=exp(j*phi);
    
end
plot(sumrate)



