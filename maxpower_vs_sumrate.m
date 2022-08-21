%emin=2   n=40
figure
sr_both=[4.02 4.79 5.27 5.63 5.92];
sr_theta=[2.06 2.09 2.061 2.07 2.06];
sr_power=[3.5804 3.72 3.93 3.97 4.03];
ga=[3.82 3.9172 4.07 4.14 4.30];
no_ris=[2.11693 2.12 2.128 2.17693 2.21];
pmax=[10 20 30 40 50];
hold on
plot(pmax, sr_both,'--or','linewidth',2);
plot(pmax, ga,'--ok','linewidth',2)
plot(pmax, sr_power,'--ob','linewidth',2)
plot(pmax,no_ris,'--om','LineWidth',2)
plot(pmax, sr_theta,'--og','linewidth',2)


grid on
xticks(pmax)
xticklabels({'10','20','30','40','50'})
xlabel('P_{max} (mW)')
ylabel('Sum Rate (bps/Hz)')
legend('Proposed Algorithm','Genetic Algorithm','Maximum Power','No RIS','Random Phase','location','best')