load('comparisons_1')
hold on    
plot((final__)/idx,'r',"linewidth",1.5)
plot((final__power)/idx,'b',"linewidth",1.5)
plot((final__theta)/idx,'g',"linewidth",1.5)
sum(x_)/length(x_)
xlabel("Iterations")
ylabel("Sum Rate")
grid on
legend('Theta+Power Optimized','Maximum Power','Random Phase','location','best')