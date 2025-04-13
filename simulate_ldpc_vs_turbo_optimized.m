%% simulate_ldpc_vs_turbo_optimized.m
% 本脚本对比 LDPC 与 Turbo 码在 AWGN 信道下的性能（n = 648, rate = 0.5）
% 固定迭代次数为 50 次，并增加了测量两种码的解码时间。
%
% LDPC 部分：
%   - 使用系统型奇偶校验矩阵 H = [H1 | I]，其中 H1 每列包含 4 个 1
%   - 译码时固定进行 50 次迭代
%
% Turbo 部分：
%   - 使用 MATLAB 内置的 TurboEncoder/Decoder 对象，采用默认 Trellis (poly2trellis(4, [13 15], 13))
%   - 固定迭代次数为 50 次，采用 'Max*' 算法，并设定 puncturing 参数
%
% 仿真中采用 BPSK 调制、AWGN 信道，通过 BER-SNR 曲线展示性能和平均解码时间。

%% 仿真参数
snr_dB    = 0:1:10;     % SNR 范围 (dB)
numFrames = 200;        % 每个 SNR 点的仿真帧数
k         = 324;        % 信息位数
N_ldpc    = 648;        % LDPC 码字长度 (码率 = 324/648 = 0.5)

%% LDPC 码构造
M_ldpc = N_ldpc - k;    % 奇偶校验矩阵的行数
d_v = 4;                % 每列 1 的个数设为 4（增强连通性）
H1 = generateH1(M_ldpc, M_ldpc, d_v);  % 生成稀疏 H1 矩阵
I_part = sparse(logical(speye(M_ldpc))); % 单位矩阵
H_ldpc = [H1, I_part];  % 系统型奇偶校验矩阵

% 创建 LDPC 编码器与译码器配置对象
encoderConfig_ldpc = ldpcEncoderConfig(H_ldpc);
decoderConfig_ldpc = ldpcDecoderConfig(encoderConfig_ldpc);
maxLDPCIter = 50;  % 固定译码迭代次数

%% Turbo 码构造
% 使用默认 Trellis: poly2trellis(4, [13 15], 13)
trellisTurbo = poly2trellis(4, [13 15], 13);
pMLen = log2(trellisTurbo.numStates);      % 尾部比特数，通常为3
pN = log2(trellisTurbo.numOutputSymbols);    % 每个编码器输出的比特数，通常为2
fullOutputLen = (k + pMLen) * pN * 2;         % Turbo 编码全输出长度，例如(324+3)*2*2 = 1308

% 设计 puncturing 索引：均匀选取 648 个符号实现率 1/2
tempIdx = (1:2:1296).';    % 生成 648 个索引
outIndices = tempIdx;      % 用于 TurboEncoder 的 OutputIndices
% TurboDecoder中使用相同的 puncturing pattern，通过 InputIndices 指定
% 设定交织器：这里采用逆序排列
interleaverIndices = (k:-1:1).';

% 创建 Turbo 编码器
turboEnc = comm.TurboEncoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'OutputIndicesSource', 'Property', ...
    'OutputIndices', outIndices);

% 创建 Turbo 译码器：固定 50 次迭代，使用 'Max*' 算法及 puncturing 信息
turboDec = comm.TurboDecoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'NumIterations', 50, ...
    'Algorithm', 'Max*', ...
    'NumScalingBits', 3, ...
    'InputIndicesSource', 'Property', ...
    'InputIndices', outIndices);

%% 定义 BPSK 调制函数
% 映射规则: 0 -> +1, 1 -> -1
bpskMod = @(bits) 1 - 2*double(bits);

%% 初始化 BER 与解码时间统计变量
ber_ldpc   = zeros(size(snr_dB));
ber_turbo  = zeros(size(snr_dB));
decTime_ldpc = zeros(size(snr_dB));   % 平均 LDPC 译码时间（秒）
decTime_turbo = zeros(size(snr_dB));  % 平均 Turbo 译码时间（秒）

%% 仿真循环
fprintf('开始仿真（优化后）：LDPC 与 Turbo 码对比（含解码时间测量）\n');
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    numErrors_ldpc = 0;  numTotal_ldpc = 0;
    numErrors_turbo = 0; numTotal_turbo = 0;
    totalDecTime_ldpc = 0;  % 累计 LDPC 译码时间
    totalDecTime_turbo = 0; % 累计 Turbo 译码时间
    
    % 对于 BPSK, 噪声方差 noiseVar = 1/(2*10^(snr/10))
    noiseVar = 1/(2*10^(snr/10));
    
    for frame = 1:numFrames
        %% 生成随机信息比特 (324x1)
        infoBits = randi([0 1], k, 1) > 0;
        
        %% LDPC 编码、调制及传输（AWGN信道）
        codeword_ldpc = ldpcEncode(infoBits, encoderConfig_ldpc);
        txSymbols_ldpc = bpskMod(codeword_ldpc);
        rxSymbols_ldpc = txSymbols_ldpc + sqrt(noiseVar)*randn(N_ldpc,1);
        rxLLR_ldpc = 2*rxSymbols_ldpc ./ noiseVar;
        
        % 测量 LDPC 译码时间
        tStart = tic;
        decodedBits_ldpc = ldpcDecode(rxLLR_ldpc, decoderConfig_ldpc, maxLDPCIter, ...
                          'OutputFormat', 'info', 'DecisionType', 'hard');
        tElapsed = toc(tStart);
        totalDecTime_ldpc = totalDecTime_ldpc + tElapsed;
        
        numErrors_ldpc = numErrors_ldpc + sum(infoBits ~= decodedBits_ldpc);
        numTotal_ldpc = numTotal_ldpc + k;
        
        %% Turbo 编码、调制及传输（AWGN信道）
        encoded_turbo = step(turboEnc, infoBits);  % 输出长度应为 648
        txSymbols_turbo = bpskMod(encoded_turbo);
        rxSymbols_turbo = txSymbols_turbo + sqrt(noiseVar)*randn(length(encoded_turbo),1);
        rxLLR_turbo = (-2/(noiseVar/2)) * real(rxSymbols_turbo);
        
        % 测量 Turbo 译码时间
        tStart = tic;
        decodedBits_turbo = step(turboDec, rxLLR_turbo);
        tElapsed = toc(tStart);
        totalDecTime_turbo = totalDecTime_turbo + tElapsed;
        
        numErrors_turbo = numErrors_turbo + sum(infoBits ~= decodedBits_turbo);
        numTotal_turbo = numTotal_turbo + k;
    end
    
    ber_ldpc(s)   = numErrors_ldpc / numTotal_ldpc;
    ber_turbo(s)  = numErrors_turbo / numTotal_turbo;
    
    % 计算平均译码时间（单位：秒）
    decTime_ldpc(s) = totalDecTime_ldpc / numFrames;
    decTime_turbo(s) = totalDecTime_turbo / numFrames;
    
    fprintf('SNR = %.1f dB: LDPC BER = %e (Avg Dec Time = %.4f s), Turbo BER = %e (Avg Dec Time = %.4f s)\n', ...
        snr, ber_ldpc(s), decTime_ldpc(s), ber_turbo(s), decTime_turbo(s));
end

%% 绘制 BER-SNR 比较图
figure;
semilogy(snr_dB, ber_ldpc, 'b-o', 'LineWidth',1.5);
hold on;
semilogy(snr_dB, ber_turbo, 'r-s', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('LDPC vs Turbo 码 (n = 648, rate = 0.5) 在 AWGN 信道下的性能与平均译码时间');
legend('LDPC', 'Turbo');
hold off;

%% 可选：绘制平均译码时间对比图
figure;
plot(snr_dB, decTime_ldpc, 'b-o', 'LineWidth',1.5);
hold on;
plot(snr_dB, decTime_turbo, 'r-s', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Average Decoding Time (s)');
title('平均译码时间: LDPC vs Turbo');
legend('LDPC', 'Turbo');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 辅助函数：生成 H1 矩阵
function H1 = generateH1(numRows, numCols, d_v)
% generateH1 生成大小为 (numRows x numCols) 的稀疏逻辑型矩阵，
% 每列恰好含有 d_v 个 1，其余为 0.
H1 = false(numRows, numCols);
for j = 1:numCols
    rows = randperm(numRows, d_v);
    H1(rows, j) = true;
end
H1 = sparse(H1);
end