%% simulate_ldpc_vs_turbo_optimized_Rayleigh.m
% 本脚本对比 LDPC 与 Turbo 码在瑞利衰落信道下的性能，二者均采用：
%   码长 n = 648, 信息位数 k = 324, 码率 = 0.5, 固定迭代次数均为 50 次。
%
% 主要修改：
%   1. 信道从 AWGN 信道改为瑞利衰落信道，发送信号先乘以瑞利衰落系数，
%      再加上复数高斯噪声；接收端进行相干检测（假设完美信道估计）。
%   2. 对 LDPC 部分保持 H1 每列 4 个 1，Turbo 码采用 'Max*' 算法并固定 50 次迭代，
%      同时 puncturing 策略保持不变。
%
% 两种编码均采用 BPSK 调制，仿真后绘制出 BER-SNR 曲线。

%% 仿真参数设置
snr_dB    = 0:1:10;     % SNR 范围，单位 dB
numFrames = 200;        % 每个 SNR 点的仿真帧数（可根据需要调整）
k         = 324;        % 信息位数
N_ldpc    = 648;        % LDPC 码字长度（率 = 324/648 = 0.5）

%% LDPC 码构造
% 构造 LDPC 奇偶校验矩阵 H = [H1 | I]
M_ldpc = N_ldpc - k;    % 奇偶校验矩阵行数
d_v = 4;                % 每列 1 的个数设为 4（优化提取连接度）
H1 = generateH1(M_ldpc, M_ldpc, d_v);  % 生成稀疏 H1 矩阵
I_part = sparse(logical(speye(M_ldpc))); % 单位矩阵
H_ldpc = [H1, I_part];  % 合成系统型奇偶校验矩阵

% 创建 LDPC 编码器与译码器配置对象
encoderConfig_ldpc = ldpcEncoderConfig(H_ldpc);
decoderConfig_ldpc = ldpcDecoderConfig(encoderConfig_ldpc);
maxLDPCIter = 50;       % 固定 50 次迭代

%% Turbo 码构造
% 使用默认 Trellis：poly2trellis(4, [13 15], 13)
trellisTurbo = poly2trellis(4, [13 15], 13);
pMLen = log2(trellisTurbo.numStates);      % 通常为 3
pN = log2(trellisTurbo.numOutputSymbols);    % 通常为 2
fullOutputLen = (k + pMLen) * pN * 2;         % 例如 (324+3)*2*2 = 1308

% 设计 puncturing 索引：均匀选取 648 个符号实现率 1/2
tempIdx = (1:2:1296).';    % 生成 648 个索引
outIndices = tempIdx;      % 用作 TurboEncoder 的 OutputIndices
% 对 TurboDecoder，puncturing 信息通过 InputIndices 指定，同样使用 outIndices

% 设定简单的交织器（逆序排列作为示例）
interleaverIndices = (k:-1:1).';

% 创建 Turbo 编码器
turboEnc = comm.TurboEncoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'OutputIndicesSource', 'Property', ...
    'OutputIndices', outIndices);

% 创建 Turbo 译码器：固定 50 次迭代，采用 'Max*' 算法及适当缩放，
% 并将 puncturing 信息通过 InputIndices 指定。
turboDec = comm.TurboDecoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'NumIterations', 50, ...
    'Algorithm', 'Max*', ...
    'NumScalingBits', 3, ...
    'InputIndicesSource', 'Property', ...
    'InputIndices', outIndices);

%% 定义 BPSK 调制映射函数
% 映射规则： 0 -> +1, 1 -> -1
bpskMod = @(bits) 1 - 2*double(bits);

%% 初始化 BER 储存变量
ber_ldpc   = zeros(size(snr_dB));
ber_turbo  = zeros(size(snr_dB));

%% 仿真循环
fprintf('开始仿真（瑞利衰落信道下）：LDPC 与 Turbo 码对比\n');
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    numErrors_ldpc = 0;  numTotal_ldpc = 0;
    numErrors_turbo = 0; numTotal_turbo = 0;
    
    % 对于 BPSK，相对于 AWGN 信道，此处 SNR 表示每比特信噪比，
    % 噪声方差 noiseVar = 1/(2*10^(snr/10))
    noiseVar = 1/(2*10^(snr/10));
    
    for frame = 1:numFrames
        %% 随机生成信息比特（长度为 324）
        infoBits = randi([0 1], k, 1) > 0;
        
        %% LDPC 编码、调制及瑞利衰落信道传输
        % LDPC 编码
        codeword_ldpc = ldpcEncode(infoBits, encoderConfig_ldpc);
        % BPSK 调制
        txSymbols_ldpc = bpskMod(codeword_ldpc);
        
        % 瑞利衰落系数：复数高斯（归一化）
        fade_ldpc = sqrt(0.5)*(randn(N_ldpc,1) + 1i*randn(N_ldpc,1));
        % 生成复数 AWGN 噪声，噪声方差按每个维度 noiseVar/2
        noise_ldpc = sqrt(noiseVar/2)*(randn(N_ldpc,1) + 1i*randn(N_ldpc,1));
        % 瑞利信道传输：乘以 fading，再加噪声
        rxSymbols_ldpc = fade_ldpc .* txSymbols_ldpc + noise_ldpc;
        % 相干检测：假设完美信道估计，对接收信号进行均衡
        rxEqualized_ldpc = rxSymbols_ldpc ./ fade_ldpc;
        % 计算 LLR（由于 txSymbols 为实数，取实部计算）
        rxLLR_ldpc = 2*real(rxEqualized_ldpc) ./ noiseVar;
        
        % LDPC 译码
        decodedBits_ldpc = ldpcDecode(rxLLR_ldpc, decoderConfig_ldpc, maxLDPCIter, ...
                          'OutputFormat', 'info', 'DecisionType', 'hard');
        numErrors_ldpc = numErrors_ldpc + sum(infoBits ~= decodedBits_ldpc);
        numTotal_ldpc = numTotal_ldpc + k;
        
        %% Turbo 编码、调制及瑞利衰落信道传输
        % Turbo 编码
        encoded_turbo = step(turboEnc, infoBits);  % 输出长度应为 648
        % BPSK 调制
        txSymbols_turbo = bpskMod(encoded_turbo);
        
        % 瑞利衰落：生成 fading 系数和噪声（长度与 encoded_turbo 相同）
        fade_turbo = sqrt(0.5)*(randn(length(encoded_turbo),1) + 1i*randn(length(encoded_turbo),1));
        noise_turbo = sqrt(noiseVar/2)*(randn(length(encoded_turbo),1) + 1i*randn(length(encoded_turbo),1));
        rxSymbols_turbo = fade_turbo .* txSymbols_turbo + noise_turbo;
        rxEqualized_turbo = rxSymbols_turbo ./ fade_turbo;
        % 计算 Turbo 部分的 LLR（同样取实部）
        rxLLR_turbo = (-2/(noiseVar/2)) * real(rxEqualized_turbo);
        
        % Turbo 译码
        decodedBits_turbo = step(turboDec, rxLLR_turbo);
        numErrors_turbo = numErrors_turbo + sum(infoBits ~= decodedBits_turbo);
        numTotal_turbo = numTotal_turbo + k;
    end
    
    ber_ldpc(s)   = numErrors_ldpc / numTotal_ldpc;
    ber_turbo(s)  = numErrors_turbo / numTotal_turbo;
    
    fprintf('SNR = %.1f dB: LDPC BER = %e, Turbo BER = %e\n', snr, ber_ldpc(s), ber_turbo(s));
end

%% 绘制 BER-SNR 性能对比图
figure;
semilogy(snr_dB, ber_ldpc, 'b-o', 'LineWidth',1.5);
hold on;
semilogy(snr_dB, ber_turbo, 'r-s', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('LDPC vs Turbo 码 (n = 648, rate = 0.5) 在瑞利衰落信道下的性能');
legend('LDPC', 'Turbo');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 辅助函数：生成 H1 矩阵
function H1 = generateH1(numRows, numCols, d_v)
% generateH1 生成大小为 (numRows x numCols) 的稀疏逻辑型矩阵，
% 每列包含 d_v 个 1，其余位置为 0.
H1 = false(numRows, numCols);
for j = 1:numCols
    rows = randperm(numRows, d_v);
    H1(rows, j) = true;
end
H1 = sparse(H1);
end