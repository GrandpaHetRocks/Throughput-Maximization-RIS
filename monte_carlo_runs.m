clc
clear all
% rng(325)

final__=zeros(1,51);

run=1;

%  k=4;
%     n=10;
%     P=ones(1,k)*5;
%     a=randn(k,k);
%     b=randn(k,k);
%     hij_array=complex(a,b);
%     a=randn(1,n);
%     b=randn(1,n);
%     theta=complex(a,b)';
%     clear sum
%     theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));
%     Hri_array=[];
%     for i=1:k
%         a=randn(n,1);
%         b=randn(n,1);
%         Hri=complex(a,b);
%         Hri_array=[Hri_array Hri];
%     end
% 
%     a=randn(1,k);
%     b=randn(1,k);
%     G_array=complex(a,b);


scaling=10^6;

for runs=1:run
    k=4;
    n=100;
    P=ones(1,k)/scaling;
    [hij_array,theta,Hri_array,G_array]=channel_generation(k,n);
    alphak=0;
%     sumrate1=[];
%     sumrate2=[];
%     sumrate3=[];
%     sumrate4=[];

    betas=zeros(k,1);
    alphas=zeros(k,1);
%     Us=zeros(k,1);
        U1=complex(zeros(n,n),zeros(n,n));
        U2=complex(zeros(n,n),zeros(n,n));
        U3=complex(zeros(n,n),zeros(n,n));
        U4=complex(zeros(n,n),zeros(n,n));
    vs=[complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1))];

    grad_=[];
    final=[];
%     for i=1:k
%             sum3=0;
%             for j=1:k
%                 if(j~=i)
%                     sum3=sum3+abs(hij_array(i,j)'+G_array(j)*theta'*Hri_array(:,j))^2*P(j);
%                 end
%             end
%             gamma=abs(hij_array(i,i)'+G_array(i)*theta'*Hri_array(:,i))^2 * P(i)/sum3;
%             if(i==1)
%             sumrate1=[sumrate1 log2(1+gamma)];
%             end
%             if(i==2)
%             sumrate2=[sumrate2 log2(1+gamma)];
%             end
%             if(i==3)
%             sumrate3=[sumrate3 log2(1+gamma)];
%             end
%             if(i==4)
%             sumrate4=[sumrate4 log2(1+gamma)];
%             end
% 
%     end


%     net_sumrate=(sumrate1+sumrate2+sumrate3+sumrate4);
    net_sumrate=sumrate_calc(hij_array,theta,Hri_array,G_array,k,P);%,sumrate1,sumrate2,sumrate3,sumrate4);
    final=[final net_sumrate];

        for iterates=1:50
            for i=1:k
                sum4=0;
                for j=1:k
%                     if(j~=i)
                    hk_=P(j)*abs(hij_array(i,j)'+G_array(j)*theta'*Hri_array(:,i))^2;
                    sum4=sum4+hk_;
%                     end
                end
                betai=(P(i)*sqrt((1+alphas(i))) * (hij_array(i,i)'+ G_array(i)*theta'*Hri_array(:,i)))/sum4;
                betas(i)=betai;
                hi=hij_array(i,i)+G_array(i)*Hri_array(:,i)'*theta;
                kaii=1/(P(i))* real(conj(betas(i))*hi'*P(i));
                alphai=(kaii^2+kaii*sqrt(kaii^2+4))/2;
                alphas(i)=alphai;
            end

            % work on U and v


            for i=1:k
                sum4=0;
                for j=1:k
                    aa=G_array(j)*Hri_array(:,i)*P(j);
                    sum4=sum4+aa*aa';  %%%%%%CHANGE%%%%%%
%                     Us(i)=sum;
                    if(i==1)
                        U1=sum4;
                    end
                    if(i==2)
                        U2=sum4;
                    end
                    if(i==3)
                        U3=sum4;
                    end
                    if(i==4)
                        U4=sum4;
                    end
                end
            end

            U=0;
%             for i=1:k
%                 U=U+(abs(betas(i))^2)*Us(i);
%             end

               U=U1+U2+U3+U4;
%                U=U1;

%             U=U/k;

            
            for i=1:k
                sum4=0;
                for j=1:k
                    aa=G_array(j)*Hri_array(:,i)*P(j);
                    bb=conj(hij_array(i,j))*P(j);
                    sum4=sum4+ bb*aa;
                    vs(:,i)=sum4;
                end
            end

            v=0;
            for i=1:k
                aii=G_array(i)*Hri_array(:,i)*P(i);
                v=v+ (P(i)*sqrt(1+alphas(i))*conj(betas(i))*aii) - (abs(betas(i))^2)*vs(i);
            end

%             v=v/k;

            phi=angle(theta);
            
            oof=U*exp(j.*phi)-v;
            oof2=-j.*exp(-j.*phi);
            grad=2*real(oof.*oof2);
            grad_=[grad_ grad];

        %     yes=1;
%             for stepp=1:10:10000
%                 step=1/(stepp);
%                 phi1=phi-grad/step;
%                 f1=exp(j.*phi)'*U*exp(j.*phi)-2*real(v'*exp(j.*phi));
%                 f2=exp(j.*phi1)'*U*exp(j.*phi1)-2*real(v'*exp(j.*phi1));
%                 if(f1-f2>=((grad'*grad)))
% % %                     stepp
%                     break;
%                 end
%             end
              diagtheta=diag(theta');
              Ltheta=SCA_phi_step_para(U,v,n,diagtheta)*scaling;
%               Ltheta=length(theta);
              rho0=1/(Ltheta)*100;
              rhom=0.95;
              m=0;
              sig=0.4;
              while(1)
    %             step=1/(start);
                rho=rho0*rhom^m;
                phi1=phi-grad*rho;
                f1=exp(j.*phi)'*U*exp(j.*phi)-2*real(v'*exp(j.*phi));
                f2=exp(j.*phi1)'*U*exp(j.*phi1)-2*real(v'*exp(j.*phi1));
                if(f1-f2>(-grad'*grad)*sig*rho)
                    break
                end
                if(rho<1/Ltheta/10)
                    break
                end
                m=m+1;
    %             start=start-0.2*start;
              end
    %         step=100000;
%             phi2=phi-grad*step;
            phi2=phi-grad*rho;
%                     grad
%                 phi2=phi-grad;
            theta2=exp(j.*phi2);

            for i=1:k
                sum3=0;
                for j=1:k
                    if(j~=i)
                        sum3=sum3+abs(hij_array(i,j)'+G_array(j)*theta2'*Hri_array(:,j))^2*P(j);
                    end
                end
                gamma=abs(hij_array(i,i)'+G_array(i)*theta2'*Hri_array(:,i))^2 * P(i)/sum3;
                if(i==1)
                sumrate1=log2(1+gamma);
                end
                if(i==2)
                sumrate2=log2(1+gamma);
                end
                if(i==3)
                sumrate3=log2(1+gamma);
                end
                if(i==4)
                sumrate4=log2(1+gamma);
                end

            end

            new_net=(sumrate1+sumrate2+sumrate3+sumrate4);

            if(new_net>net_sumrate)
                net_sumrate=new_net;
            end
%             net_sumrate=new_net;
            phi=phi2;
            theta=exp(j.*phi);
            final=[final net_sumrate];   
        end
    final__= final__ + final;
end

hold on    
plot((final__)/run,"linewidth",1.5)
xlabel("Iterations")
ylabel("Sum Rate")

[x,fval]=ga(@(net_sumrate) sumrate_calc(hij_array,theta,Hri_array,G_array,k,P),7)




