sr_both=[6.27527 5.85672 4.91803 3.58758 3.42081 3.009 2.67902];
sr_theta=[2.23 2.20 2.18 2.27 2.19 2.25 2.25];
sr_power=[6.0579 5.30 4.18963 3.0279 2.78823 2.581 2.25086];
ga=[5.2364 4.7032 3.9210 3.0228 2.7882 2.2509 2.25];
emin=[100 500 1000 1500 3000 4000 5000];
hold on
plot(emin, sr_both,'--or','linewidth',2);
plot(emin, sr_power,'--ob','linewidth',2)
plot(emin, ga,'--ok','linewidth',2)
plot(emin, sr_theta,'--og','linewidth',2)

grid on
xticks(emin)
xticklabels({'100','500','1000','1500','3000','4000','5000'})
xlabel('emin')
ylabel('Sum Rate')
legend('Proposed Algorithm','Maximum Power','Genetic Algorithm','Random Phase','location','best')