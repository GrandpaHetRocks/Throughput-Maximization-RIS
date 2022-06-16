sr_both=[2.80064 4.24896 5.47434 5.88465 6.11297 6.2757];
sr_theta=[2.23 2.23485 2.23485 2.23485 2.23 2.23485];
sr_power=[2.20903 4.05508 5.2569 5.66722 5.89553 6.0578];
pmax=[1 10 50 100 150 200];
hold on
plot(pmax, sr_both,'--or','linewidth',2);
plot(pmax, sr_theta,'--og','linewidth',2)
plot(pmax, sr_power,'--ob','linewidth',2)
grid on
xticks(pmax)
xticklabels({'1','10','50','100','150','200'})
xlabel('pmax')
ylabel('Sum Rate')
legend('Proposed Algorithm','Random Phase','Maximum Power')