%peripheral_agent
x=8:0.01:20;
curve1=(1-x+sqrt(x.^2-10*x+17))/2;
curve2=(1-x-sqrt(x.^2-10*x+17))/2;
plot(x, curve1, 'r', 'LineWidth', 2);
hold on;
plot(x, curve2, 'b', 'LineWidth', 2);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'o');