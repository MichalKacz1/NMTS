a = 0;
b = 3;
A = 0.02;
% Przypadek a
z0_a = -abs(b - a);
G_a = tf([1, -z0_a], conv([1 a], [1 b]));

z0_b = 0.5 * (a + b);
G_b = tf([1, -z0_b], conv([1 a], [1 b]));

figure;
bode(G_a);
title('Charakterystyka Bodego - Przypadek a: z_0 = -|b - a|');
ylabel('Amplituda (dB)');
xlabel('Częstotliwość (rad/s)');
grid on;

figure;
bode(G_b);
title('Charakterystyka Bodego - Przypadek b: z_0 = 0.5(b + a)');
ylabel('Amplituda (dB)');
xlabel('Częstotliwość (rad/s)');
grid on;