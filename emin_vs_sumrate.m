%power=50 n=40
figure
sr_both=[5.92 5.07 4.43 3.59 2.84];
sr_theta=[2.06 2.17 2.04 1.94 1.79];
sr_power=[4.03 3.33 3.06 2.77 2.0637];
ga=[4.30 3.60 3.504 2.87 2.2392];
no_ris=[2.21 2.23 2.20 2.18 2.16];
emin=[2 3 4 5 6]*10^-6;
hold on
plot(emin, sr_both,'--or','linewidth',2);
plot(emin, ga,'--ok','linewidth',2)
plot(emin, sr_power,'--ob','linewidth',2)
plot(emin,no_ris,'--om','LineWidth',2)
plot(emin, sr_theta,'--og','linewidth',2)

grid on
xticks(emin)
% xticklabels({'2','3','4','5','6'})
xlabel('E_{min}')
ylabel('Sum Rate (bps/Hz)')
legend('Proposed Algorithm','Genetic Algorithm','Maximum Power','No RIS','Random Phase','location','best')