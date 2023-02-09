popv = populations1(1,:);
popw = populations1(2,:);
time = linspace(0,1000,2000);

figure(1)
plot(time, popv, LineWidth=2)
hold on
plot(time,popw, LineWidth=2)
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
title('The population rates', 'Interpreter', 'latex')
xlabel('Time','Interpreter','latex');ylabel('Population rates','Interpreter','latex');
legend('Population rate for preys', 'Population rate for predators', 'Interpreter','latex','Location','best');
grid on
axis padded
figure(2)
plot(popv,popw, LineWidth=2)
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
title('The phase portrait', 'Interpreter', 'latex')
xlabel('Population rate of preys','Interpreter','latex');ylabel('Population rates of predators','Interpreter','latex');
grid on
axis padded