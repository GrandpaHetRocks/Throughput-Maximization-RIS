% load('comparisons_1.mat')
hold on    
plot([0:50],(final__)/idx,'r',"linewidth",1.5)
plot([0:50],(final__power)/idx,'b',"linewidth",1.5)
plot([0:50],(final_whoo)/idx_,'m',"linewidth",1.5)
plot([0:50],(final__theta)/idx,'g',"linewidth",1.5)
sum(x_)/length(x_)
xlabel("Iterations")
ylabel("Sum Rate (bps/Hz)")
grid on
legend('Proposed Algorithm','Maximum Power','No RIS','Random Phase','location','best')