clear all; close all; clc;

a = 5.43e-10;
m0 = 9.1e-31;
h = 6.63e-34;

V = @(x)1e-19 * cos(2 * pi / a .* x) .* (x >= -a / 4 & x <= a / 4);

x = linspace(-a / 2, a / 2, 1000);
plot(x, V(x));
title('V(x)');
xlabel('x/m');
ylabel('V/J');

N = 8;  % for simplicity, let N be an even integer.
assert(mod(N, 2) == 0);
VN = zeros(1, 2 * N + 1);
for i = 1: 1: 2 * N + 1
    n = i - N - 1;
    VN(i) = real(integral(@(x)V(x) .* exp(-1j * n .* x * 2 * pi / a), -a / 2, a / 2) / a);
end
figure(2);
stem(-N: N, VN);
xlim([-N, N]);
title('Fourier Series');
xlabel('N');

VN_mat = zeros(N + 1);
for i = 1: 1: N + 1
    VN_mat(i, :) = VN(N + i: -1: i);
end
cof = (h / 2 / pi) ^ 2 / 2 / m0;
sample_num = 1001;  % must be odd
assert(mod(sample_num, 2) == 1);
k0 = linspace(-pi / a, pi / a, sample_num);
eigen_E = zeros(sample_num, N + 1);
for i = 1: 1: sample_num
    kN_val = linspace(k0(i) - 2 * pi * N / 2 / a, k0(i) + 2 * pi * N / 2 / a, N + 1);
    offset_diag = diag(kN_val .^ 2 .* cof);
    eig_mat = VN_mat + offset_diag;
    eigen_E(i, :) = eig(eig_mat);
end

figure(3);
draw_num = N + 1;
assert(draw_num <= N + 1);
for i = 1: 1: draw_num
    plot(k0, eigen_E(:, i));
    xlim([k0(1), k0(end)]);
    if i ~= draw_num
        hold on;
    else
        title('Reduced Brillouin Zone');
        xlabel('k');
        ylabel('E');
    end
end

real_energy_gap = min(eigen_E(:, 2: end) - eigen_E(:, 1: end - 1), [], 1);
est_energy_gap = 2 * abs(VN(N + 2: 2 * N + 1));
figure(4)
scatter(1: 1: N, real_energy_gap);
hold on
scatter(1: 1: N, est_energy_gap);
xlim([0, N]);
legend('Real', 'Est');
title('Energy Gap');
xlabel('N');
ylabel('E');