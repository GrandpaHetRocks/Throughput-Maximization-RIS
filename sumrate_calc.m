function [net_sumrate]=sumrate_calc(hij_array,theta,Hri_array,G_array,k,P,sumrate1,sumrate2,sumrate3,sumrate4)
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

    net_sumrate=norm((sumrate1+sumrate2+sumrate3+sumrate4));
end