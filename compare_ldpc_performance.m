%% compare_ldpc_performance.m
% 本脚本加载 simulate_short_ldpc_data.m 保存的 results.mat 数据，
% 并在同一图中对比不同最大译码迭代次数下的 BER-SNR 性能。

clear; clc;
load('results.mat', 'results');

snr_dB = results.snr_dB;
iterCounts = results.iterCounts;
numIterations = length(iterCounts);

figure;
hold on;
markers = {'o', 's', 'd', '^', 'v'};  % 不同曲线使用不同的marker
for i = 1:numIterations
    ber_curve = results.data{i};
    plot(snr_dB, ber_curve, ['-b' markers{mod(i-1, length(markers))+1}], 'LineWidth',1.5);
end
grid on;
xlabel('SNR (dB)');
ylabel('比特错误率 (BER)');
title('不同译码最大迭代次数下的 LDPC 码性能对比');
legendStrings = cell(numIterations,1);
for i = 1:numIterations
    legendStrings{i} = sprintf('maxIter = %d', iterCounts(i));
end
legend(legendStrings, 'Location', 'northeast');
hold off;