%power=50   emin=4
figure
sr_both=[3.58 3.7194 5.52 5.61];
sr_theta=[1.8769 1.90 2.09 2.26];
sr_power=[2.66 2.89 4.04 4.45];
ga=[2.68 3.40 4.53 5.28];
no_ris=[2.24 2.24 2.24 2.24];
n=[10 30 50 70];
hold on
plot(n, sr_both,'--or','linewidth',2);
plot(n, ga,'--ok','linewidth',2)
plot(n, sr_power,'--ob','linewidth',2)
plot(n, no_ris,'--om','LineWidth',2)
plot(n, sr_theta,'--og','linewidth',2)

grid on
xticks(n)
xticklabels({'10','30','50','70'})
xlabel('N')
ylabel('Sum Rate (bps/Hz)')
legend('Proposed Algorithm','Genetic Algorithm','Maximum Power','No RIS','Random Phase','location','best')