%% simulate_ldpc_vs_turbo_optimized.m
% 本脚本对比 LDPC 与 Turbo 码在相同码长 (n = 648) 和码率 (0.5) 下的性能。
% 优化目标：
%   1. 固定两种编码的迭代次数均为 50 次。
%   2. 对 Turbo 码进行参数优化，采用 "Max*" 算法及适当 puncturing，
%      使得在低 SNR (0~1 dB) 时 Turbo 的误码率不明显高于 LDPC，
%      同时在中高 SNR 下 LDPC 码表现出误码率快速下降、接近香农极限的特性。
%
% 两种编码均采用 BPSK 调制和 AWGN 信道，并统计 BER 曲线。

%% 仿真参数设置
snr_dB    = 0:1:10;     % SNR (dB) 范围，从 0 dB 到 10 dB
numFrames = 200;        % 每个 SNR 下的仿真帧数（可根据需要增加以提高统计精度）
k         = 324;        % 信息位数
N_ldpc    = 648;        % LDPC 码字长度 (码率 = 324/648 = 0.5)

%% LDPC 码构造
% 1. 构造 LDPC 奇偶校验矩阵 H = [H1 | I]
% 为增强 LDPC 码的性能，将 H1 部分中每列 1 的个数由 3 增加到 4
M_ldpc = N_ldpc - k;    % 奇偶校验矩阵的行数
d_v = 4;                % 每列 1 的个数设为 4
H1 = generateH1(M_ldpc, M_ldpc, d_v);  % 生成稀疏 H1 矩阵
I_part = sparse(logical(speye(M_ldpc))); % 单位矩阵
H_ldpc = [H1, I_part];  % 合成系统型奇偶校验矩阵

% 2. 构造 LDPC 编码器与译码器
encoderConfig_ldpc = ldpcEncoderConfig(H_ldpc);
decoderConfig_ldpc = ldpcDecoderConfig(encoderConfig_ldpc);
maxLDPCIter = 50;       % 固定 LDPC 译码迭代次数为 50

%% Turbo 码构造
% 使用 MATLAB 内置的 TurboEncoder 和 TurboDecoder 对象，并对参数进行优化：
% - 采用默认 Trellis (poly2trellis(4, [13 15], 13))
% - 使用简单均匀的 puncturing 策略实现码率 1/2
% - 将 Turbo 译码迭代次数固定为 50，同时选择 'Max*' 算法改善低 SNR 性能
trellisTurbo = poly2trellis(4, [13 15], 13);
pMLen = log2(trellisTurbo.numStates);      % tail bits 数 = log2(numStates)，通常为 3
pN = log2(trellisTurbo.numOutputSymbols);    % 每个编码器输出比特数，通常为 2
fullOutputLen = (k + pMLen) * pN * 2;         % Turbo 编码全输出长度 (例如: (324+3)*2*2 = 1308)

% 设计 puncturing 索引：选取 fullOutputLen 内均匀的部分使输出达到 648 （码率 1/2）
tempIdx = (1:2:1296).';    % 这里产生 648 个索引
outIndices = tempIdx;      % 将此索引用于 TurboEncoder 的 OutputIndices
% 对 TurboDecoder，puncturing 信息通过 InputIndices 指定，同样使用 outIndices

% 设定简单交织器，这里采用逆序排列（仅作为示例，可根据需要设计更优的交织器）
interleaverIndices = (k:-1:1).';

% 构造 Turbo 编码器
turboEnc = comm.TurboEncoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'OutputIndicesSource', 'Property', ...
    'OutputIndices', outIndices);

% 构造 Turbo 译码器：
% 设置固定 50 次迭代，采用 'Max*' 算法和适当缩放，puncturing 信息通过 InputIndices 指定
turboDec = comm.TurboDecoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'NumIterations', 50, ...           % 固定迭代次数 50
    'Algorithm', 'Max*', ...           % 采用 'Max*' 算法增强低 SNR 性能
    'NumScalingBits', 3, ...
    'InputIndicesSource', 'Property', ...
    'InputIndices', outIndices);

%% 定义 BPSK 调制映射函数
% 映射规则： 0 -> +1, 1 -> -1
bpskMod = @(bits) 1 - 2*double(bits);

%% 初始化 BER 结果存储
ber_ldpc   = zeros(size(snr_dB));
ber_turbo  = zeros(size(snr_dB));

%% 仿真循环
fprintf('开始仿真（优化后）：LDPC 与 Turbo 码对比\n');
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    numErrors_ldpc = 0;  numTotal_ldpc = 0;
    numErrors_turbo = 0; numTotal_turbo = 0;
    
    % 对于 BPSK, 信道噪声方差 noiseVar = 1/(2*10^(snr/10))
    noiseVar = 1/(2*10^(snr/10));
    
    for frame = 1:numFrames
        %% 随机生成信息比特（长度：324）
        infoBits = randi([0 1], k, 1) > 0;
        
        %% LDPC 编码、信道传输与译码过程
        codeword_ldpc = ldpcEncode(infoBits, encoderConfig_ldpc);  % LDPC 编码
        txSymbols_ldpc = bpskMod(codeword_ldpc);                   % BPSK 调制
        rxSymbols_ldpc = txSymbols_ldpc + sqrt(noiseVar)*randn(N_ldpc,1); % AWGN 信道
        rxLLR_ldpc = 2*rxSymbols_ldpc./noiseVar;                   % 计算对数似然比 (LLR)
        decodedBits_ldpc = ldpcDecode(rxLLR_ldpc, decoderConfig_ldpc, maxLDPCIter, ...
                          'OutputFormat', 'info', 'DecisionType', 'hard'); % LDPC 译码
        numErrors_ldpc = numErrors_ldpc + sum(infoBits ~= decodedBits_ldpc);
        numTotal_ldpc = numTotal_ldpc + k;
        
        %% Turbo 编码、信道传输与译码过程
        encoded_turbo = step(turboEnc, infoBits);  % Turbo 编码，输出长度应为 648
        txSymbols_turbo = bpskMod(encoded_turbo);    % BPSK 调制
        rxSymbols_turbo = txSymbols_turbo + sqrt(noiseVar)*randn(length(encoded_turbo),1); % AWGN 信道
        % 计算 Turbo 译码器所需的 LLR 值，注意缩放因子与实现有关
        rxLLR_turbo = (-2/(noiseVar/2)) * real(rxSymbols_turbo);
        decodedBits_turbo = step(turboDec, rxLLR_turbo); % Turbo 译码
        numErrors_turbo = numErrors_turbo + sum(infoBits ~= decodedBits_turbo);
        numTotal_turbo = numTotal_turbo + k;
    end
    ber_ldpc(s)   = numErrors_ldpc / numTotal_ldpc;
    ber_turbo(s)  = numErrors_turbo / numTotal_turbo;
    
    fprintf('SNR = %.1f dB: LDPC BER = %e, Turbo BER = %e\n', snr, ber_ldpc(s), ber_turbo(s));
end

%% 绘制 BER 曲线图
figure;
semilogy(snr_dB, ber_ldpc, 'b-o', 'LineWidth',1.5);
hold on;
semilogy(snr_dB, ber_turbo, 'r-s', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('LDPC vs Turbo 码 (n = 648, rate = 0.5) 在 AWGN 信道下的性能');
legend('LDPC', 'Turbo');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 辅助函数：生成 H1 矩阵
function H1 = generateH1(numRows, numCols, d_v)
% generateH1 生成一个大小为 (numRows x numCols) 的稀疏逻辑型矩阵，
% 每列包含 d_v 个1，其他位置均为 0.
H1 = false(numRows, numCols);
for j = 1:numCols
    rows = randperm(numRows, d_v);
    H1(rows, j) = true;
end
H1 = sparse(H1);
end