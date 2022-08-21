clc
clear all
% close all
rng(12,'Twister') %12
figure

final_whoo=zeros(1,51);
% final__power=zeros(1,51);
% final__zeros(,1)=zeros(1,51);
x_=[];

run=2000;
emin=4;  
receiver_energy=[0 0 0 0];

pow=50;

errors=0;


scaling=10^6;
% maxpow=1/scaling;
maxpow=1*pow;
k=4;

idx_=0;
tic
for runs=1:run
    nn_poweronly=0;
    value1=ones(1,k)*10000;
    n=40;
    P=ones(1,k)*pow;%/scaling;
    [hij_array,theta,Hri_array,G_array,vij_array,Vri_array]=channel_generation2(k,n);
%     if(runs==12 || runs==19)
%         continue
%     end
    alphak=0;
    lagrange=ones(1,4)*1000;

    betas=zeros(k,1);
    alphas=zeros(k,1);

        U1=complex(zeros(n,n),zeros(n,n));
        U2=complex(zeros(n,n),zeros(n,n));
        U3=complex(zeros(n,n),zeros(n,n));
        U4=complex(zeros(n,n),zeros(n,n));
    vs=[complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1)) complex(zeros(n,1),zeros(n,1))];

    grad_=[];
    final=[];
    final_poweronly=[];
%     final_zeros(n,1)only=[];
    
    net_sumrate=sumrate_calc(hij_array,zeros(n,1),Hri_array,G_array,k,P);%,sumrate1,sumrate2,sumrate3,sumrate4);
    final=[final net_sumrate];
%     final_poweronly=[final_poweronly net_sumrate];
%     final_zeros(n,1)only=[final_zeros(n,1)only net_sumrate];

        for iterates=1:50
%             for i=1:k
%                 sum4=0;
%                 for j=1:k
% %                     if(j~=i)
%                     hk_=P(j)*abs(hij_array(i,j)'+G_array(:,j)'*diag(zeros(n,1))*Hri_array(:,i))^2;
%                     sum4=sum4+hk_;
% %                     end
%                 end
%                 betai=(P(i)*sqrt((1+alphas(i))) * (hij_array(i,i)'+ G_array(:,i)'*diag(zeros(n,1))*Hri_array(:,i)))/sum4;
%                 betas(i)=betai;
% %                 size(zeros(n,1)')
% %                 size(diag(Hri_array(:,i)))
% %                 size(G_array(:,i))
%                 hi=hij_array(i,i)+zeros(n,1)'*diag(Hri_array(:,i))*G_array(:,i);
%                 kaii=1/(P(i))* real(conj(betas(i))*hi'*P(i));
%                 alphai=(kaii^2+kaii*sqrt(kaii^2+4))/2;
%                 alphas(i)=alphai;
%             end

            % work on U and v


%             for i=1:k
%                 sum4=0;
%                 for j=1:k
% %                     aa=G_array(j)*Hri_array(:,i)*P(j);
%                     aa=diag(Hri_array(:,i))*G_array(:,j)*P(j);
%                     sum4=sum4+aa*aa';  %%%%%%CHANGE%%%%%%
% %                     Us(i)=sum;
%                     if(i==1)
%                         U1=sum4;
%                     end
%                     if(i==2)
%                         U2=sum4;
%                     end
%                     if(i==3)
%                         U3=sum4;
%                     end
%                     if(i==4)
%                         U4=sum4;
%                     end
%                 end
%             end
% 
%             U=0;
% %             for i=1:k
% %                 U=U+(abs(betas(i))^2)*Us(i);
% %             end
% 
%                U=U1+U2+U3+U4;
% %                U=U1;
% 
% %             U=U/k;
% 
%             
%             for i=1:k
%                 sum4=0;
%                 for j=1:k
%                     aa=G_array(j)*Hri_array(:,i)*P(j);
%                     bb=conj(hij_array(i,j))*P(j);
%                     sum4=sum4+ bb*aa;
%                     vs(:,i)=sum4;
%                 end
%             end
% 
%             v=0;
%             for i=1:k
%                 aii=G_array(i)*Hri_array(:,i)*P(i);
%                 v=v+ (P(i)*sqrt(1+alphas(i))*conj(betas(i))*aii) - (abs(betas(i))^2)*vs(i);
%             end
% 
% %             v=v/k;
% 
%             phi=angle(zeros(n,1));
%             
%             oof=U*exp(j.*phi)-v;
%             oof2=-j.*exp(-j.*phi);
%             grad=2*real(oof.*oof2);
%             grad_=[grad_ grad];
% 
%         %     yes=1;
% %             for stepp=1:10:10000
% %                 step=1/(stepp);
% %                 phi1=phi-grad/step;
% %                 f1=exp(j.*phi)'*U*exp(j.*phi)-2*real(v'*exp(j.*phi));
% %                 f2=exp(j.*phi1)'*U*exp(j.*phi1)-2*real(v'*exp(j.*phi1));
% %                 if(f1-f2>=((grad'*grad)))
% % % %                     stepp
% %                     break;
% %                 end
% %             end
%               diagzeros(n,1)=diag(zeros(n,1)');
%               Lzeros(n,1)=SCA_phi_step_para(U,v,n,diagzeros(n,1))*scaling;
% %               Lzeros(n,1)=length(zeros(n,1));
%               rho0=1/(Lzeros(n,1))*100;
%               rhom=0.95;
%               m=0;
%               sig=0.4;
%               while(1)
%     %             step=1/(start);
%                 rho=rho0*rhom^m;
%                 phi1=phi-grad*rho;
%                 f1=exp(j.*phi)'*U*exp(j.*phi)-2*real(v'*exp(j.*phi));
%                 f2=exp(j.*phi1)'*U*exp(j.*phi1)-2*real(v'*exp(j.*phi1));
%                 if(f1-f2>(-grad'*grad)*sig*rho)
%                     break
%                 end
%                 if(rho<1/Lzeros(n,1)/10)
%                     break
%                 end
%                 m=m+1;
%     %             start=start-0.2*start;
%               end
%     %         step=100000;
% %             phi2=phi-grad*step;
%             phi2=phi-grad*rho;
% %                     grad
% %                 phi2=phi-grad;
%             zeros(n,1)=exp(j.*phi2);
            
            yi=zeros(1,k);
            
            for i=1:k
                sum3=0;
                for j=1:k
                    if(j~=i)
%                         sum3=sum3+abs(hij_array(i,j)'+G_array(j)*zeros(n,1)'*Hri_array(:,j))^2*P(j);
                        sum3=sum3+abs(hij_array(i,j)')^2*P(j);
                    end
                end
                
                num1=abs(hij_array(i,i)')^2 * P(i);
                gamma=num1/sum3;
                num2=sqrt((1+gamma)*num1);
                yi(i)=num2/sum3;
                
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
            
%             sum5=zeros(1,k);
            
            for i=1:k
                sum3=0;
                sum5=0;
                for j=1:k
                    if(j~=i)
                        sum5=sum5+ (abs(hij_array(i,j)')*yi(j)^2)^2;  %no sensor
                        sum3=sum3+abs(hij_array(i,j)')^2*P(j);  %no sensor
                    end

                 

                 %energy constraint
                 
                 sum6=0;
                 sum7=0;
                 for j=1:k
                    sum6=sum6+(abs(vij_array(i,j)')^2)*scaling;
                    sum7=sum7+(abs(vij_array(i,j)')^2)*P(i)*scaling;
                 end
%                  sum7
                 sum5=sum5-lagrange(i)*sum6;
                 lagrange(i)=lagrange(i)-(sum7-emin);
                    if(lagrange(i)<0.000001)
                    lagrange(i)=0.000001;
                    end

                end
                num1=abs(hij_array(i,i)')^2 * P(i);
                gamma=num1/sum3;
%             sum5
            valuetemp=value1;
            value1(i)=(1+gamma)*abs(hij_array(i,i)')^2 * abs(yi(i))^2/(sum5);
%             maxpow
%             value1
            if((value1(i)<0.06 && value1(i)~=0))
%                 value1(i)=valuetemp(i);
                    value1(i)=0.5;
            end
%             value1
            P(i)=min(maxpow,value1(i));
            end

%             net_sumrate_poweronly=sumrate_calc(hij_array,zeros(n,1),Hri_array,G_array,k,P);
%             if(net_sumrate_poweronly>=nn_poweronly && net_sumrate_poweronly<5)
%                 nn_poweronly=net_sumrate_poweronly;
%             else
%                 net_sumrate_poweronly=nn_poweronly;
%             end
%             final_poweronly=[final_poweronly net_sumrate_poweronly];
%             
%             net_sumrate_zeros(n,1)only=sumrate_calc(hij_array,zeros(n,1),Hri_array,G_array,k,[1 1 1 1]*pow);
%             final_zeros(n,1)only=[final_zeros(n,1)only net_sumrate_zeros(n,1)only];
            
            
            new_net=(sumrate1+sumrate2+sumrate3+sumrate4);

            if(new_net>net_sumrate)
                net_sumrate=new_net;
            end
%             net_sumrate=new_net;
%             phi=phi2;
%             zeros(n,1)=exp(j.*phi);
            final=[final net_sumrate];   
        end
    if(final(35)~=final(1))
        final_whoo= final_whoo + final;
%         final__power=final__power+final_poweronly;
%         final__zeros(n,1)=final__zeros(n,1)+final_zeros(n,1)only;
%         [x,fval]=ga(@(net_sumrate) sumrate_calc(hij_array,zeros(n,1),Hri_array,G_array,k,P),7,optimoptions('ga','MaxGenerations',1000));
%         x_=[x_ fval];
        idx_=idx_+1;

        
        for i=1:k
            energy=0;
            for j=1:k
                energy=energy+(abs(vij_array(i,j)')^2)*P(i)*scaling;
       
            end
            receiver_energy(i)=receiver_energy(i)+energy;
            if(energy<emin)
                errors=errors+1;
            end
            
        end
        P
%        receiver_energy=receiver_energy+en1';  
    end
   
    runs
%     receiver_energy
end
toc/run


hold on    
plot((final_whoo)/idx_,"linewidth",1.5)
% plot((final__power)/idx,"linewidth",1.5)
% plot((final__zeros(n,1))/idx,"linewidth",1.5)
% sum(x_)/length(x_)
% xlabel("Iterations")
% ylabel("Sum Rate")
% title('Pmax=',pow)
legend('No RIS','location','southeast')

receiver_energy/idx_
errors
idx_*4

