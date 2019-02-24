loop = 10;

current = zeros(1,loop);
for i = 2:1:loop
    current(i) = a2part2c(1/i);
end
%% 
figure(11)
plot(1:1:loop,current);
title('Narrowing of Bottle-neck')
xlabel('Spacing of L*1/x')
ylabel('Current (amps)')
xlim([2 loop])