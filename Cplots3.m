popu = populations3(1,:);
popv = populations3(2,:);
popw = populations3(3,:);
time = linspace(0,1200,2400);

figure(1)
plot(time,popu, 'r')
hold on
plot(time, popv, 'g')
hold on
plot(time,popw, 'b')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
title('The population rates', 'Interpreter', 'latex')
xlabel('Time','Interpreter','latex');ylabel('Population rates','Interpreter','latex');
legend('Population rate for mutualists','Population rate for preys', 'Population rate for predators', 'Interpreter','latex','Location','best');
grid on
axis padded
figure(2)
plot3(popu,popv,popw)
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
title('The phase portrait', 'Interpreter', 'latex')
xlabel('Population rate of mutualists','Interpreter','latex');
ylabel('Population rates of preys','Interpreter','latex');
zlabel('Population rates of predators','Interpreter','latex')
grid on
axis padded
% u - mutualist, v - prey, w - predator
