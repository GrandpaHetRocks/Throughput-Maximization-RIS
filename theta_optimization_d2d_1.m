clc
clear all
rng(325)
k=4;
n=500;
P=ones(1,k)*100;
a=randn(k,k);
b=randn(k,k);
hij_array=complex(a,b);
a=randn(1,n);
b=randn(1,n);
theta=complex(a,b)';
theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));
Hri_array=[];
for i=1:k
    a=randn(n,1);
    b=randn(n,1);
    Hri=complex(a,b);
    Hri_array=[Hri_array Hri];
end

a=randn(1,k);
b=randn(1,k);
G_array=complex(a,b);

alphak=0;
sumrate1=[];
sumrate2=[];
sumrate3=[];
sumrate4=[];

betas=zeros(k,1);
alphas=zeros(k,1);
for iterates=1:50
    for i=1:k
        sum3=0;
        for j=1:k
            if(j~=i)
                sum3=sum3+abs(hij_array(i,j)'+G_array(j)*theta'*Hri_array(:,j))^2*P(j);
            end
        end
        gamma=abs(hij_array(i,i)'+G_array(i)*theta'*Hri_array(:,i))^2 * P(i)/sum3;
        if(i==1)
        sumrate1=[sumrate1 log2(1+gamma)];
        end
        if(i==2)
        sumrate2=[sumrate2 log2(1+gamma)];
        end
        if(i==3)
        sumrate3=[sumrate3 log2(1+gamma)];
        end
        if(i==4)
        sumrate4=[sumrate4 log2(1+gamma)];
        end
        
    end
%     sumrates(i)=sumrate;
    
    
    
    for i=1:k
        sum=0;
        for j=1:k
%             if(j~=i)
            hk_=abs(hij_array(i,j)'+G_array(j)*theta'*Hri_array(:,i)*P(j))^2;
            sum=sum+hk_;
%             end
        end
        betai=(sqrt(P(i)*(1+alphas(i))) * (hij_array(i,i)'+ G_array(i)*theta'*Hri_array(:,i)))/sum;
        betas(i)=betai;
        hi=hij_array(i,i)+G_array(i)*Hri_array(:,i)'*theta;
        kaii=1/sqrt(P(i))* real(conj(betas(i))*hi'*P(i));
        alphai=(kaii^2+kaii*sqrt(kaii^2+4))/2;
        alphas(i)=alphai;
    end
    
    % work on U and v
    
    sum=0;
    for i=1:k
        for j=1:k
            aa=G_array(j)*Hri_array(:,i)*P(j);
            sum=sum+aa'*aa;
        end
    end
    
    U=0;
    for i=1:k
        U=U+(abs(betas(i))^2)*sum;
    end
    
    sum=0;
    for i=1:k
        for j=1:k
            aa=G_array(j)*Hri_array(:,i)*P(j);
            bb=conj(hij_array(i,j));
            sum=sum+ bb*aa;
        end
    end
    
    v=0;
    for i=1:k
        aii=G_array(i)*Hri_array(:,i)*P(i);
        v=v+ sqrt(P(i)*(1+alphas(i))) * conj(betas(i))*aii-abs(betas(i))^2*sum;
    end
    
    phi=angle(theta);
    
    oof=U*exp(j*phi)-v;
    grad=2*real(-j*exp(-j*oof));
    
    for stepp=1:100001:5
        step=1/stepp;
        phi1=phi-grad/step;
        f1=exp(j*phi)'*U*exp(j*phi)-2*real(v'*exp(j*phi));
        f2=exp(j*phi1)'*U*exp(j*phi1)-2*real(v'*exp(j*phi1));
        if(f1-f2>=(norm(grad)^2)*step*0.2)
            break
        end
    end
%     step=100000;
    phi=phi-grad/step;
    
    theta=exp(j*phi);
            
    
end

hold on    
plot(sumrate1+sumrate2+sumrate3+sumrate4,"linewidth",1.5)
xlabel("Iterations")
ylabel("Sum Rate")




