%% simulate_ldpc_vs_turbo_optimized_burst.m
% 本脚本比较 LDPC 与 Turbo 码在突发噪声环境下的性能
% 参数设置：码长 n = 648, 信息位数 k = 324, 码率 = 0.5；固定迭代次数均为 50 次
%
% 突发噪声环境描述：
%   - 信道平时为 AWGN 信道。
%   - 每一帧有一定概率 burstProb 会出现噪声突发，
%     在突发期间，突发长度为码字长度的 burstFrac，
%     噪声方差放大 burstFactor 倍。
%
% 两种编码均采用 BPSK 调制，其仿真结果通过 BER-SNR 曲线展示。
%
% 请根据需要修改 burstProb、burstFrac 和 burstFactor 参数。

%% 仿真参数设置
snr_dB     = 0:1:10;           % SNR 范围 (dB)
numFrames  = 200;              % 每个 SNR 点的仿真帧数
k          = 324;              % 信息位数
N_ldpc     = 648;              % LDPC 码字长度 (码率 = 324/648 = 0.5)

% 突发噪声参数
burstProb   = 0.3;            % 每帧发生突发噪声的概率
burstFrac   = 0.2;            % 突发噪声覆盖码字长度的 20%
burstFactor = 10;             % 突发噪声时噪声方差放大因子

%% LDPC 码构造
% 构造系统型奇偶校验矩阵 H = [H1 | I]
M_ldpc = N_ldpc - k;          % 奇偶校验矩阵行数
d_v    = 4;                   % H1 每列 1 的个数设为 4 (增强连通性)
H1     = generateH1(M_ldpc, M_ldpc, d_v);  % 生成稀疏 H1 矩阵
I_part = sparse(logical(speye(M_ldpc)));    % 单位矩阵
H_ldpc = [H1, I_part];         % 系统型奇偶校验矩阵

% 创建 LDPC 编码器和译码器配置对象
encoderConfig_ldpc = ldpcEncoderConfig(H_ldpc);
decoderConfig_ldpc = ldpcDecoderConfig(encoderConfig_ldpc);
maxLDPCIter = 50;             % 固定 LDPC 译码迭代次数为 50 次

%% Turbo 码构造
% 使用 MATLAB 内置的 TurboEncoder/Decoder 对象
trellisTurbo = poly2trellis(4, [13 15], 13);  % 默认 Trellis 结构
pMLen      = log2(trellisTurbo.numStates);     % Tail bits 数量，通常为 3
pN         = log2(trellisTurbo.numOutputSymbols);% 每个编码器输出比特数，通常为 2
fullOutputLen = (k + pMLen) * pN * 2;           % Turbo 编码全输出长度, 如 (324+3)*2*2 = 1308

% 设计 puncturing 索引：
% 为实现率 1/2，均匀选取 648 个符号
tempIdx    = (1:2:1296).';     % 生成 648 个索引
outIndices = tempIdx;           % 对应于 TurboEncoder 的 OutputIndices
% 对 TurboDecoder，puncturing 信息通过 InputIndices 指定 (同样使用 outIndices)

% 定义简单交织器 (逆序排列作为示例)
interleaverIndices = (k:-1:1).';

% 创建 Turbo 编码器
turboEnc = comm.TurboEncoder(... 
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'OutputIndicesSource', 'Property', ...
    'OutputIndices', outIndices);

% 创建 Turbo 译码器，优化参数包括：
%   - 固定 50 次迭代，
%   - 使用 'Max*' 算法，
%   - 利用 InputIndices 指定 puncturing 信息
turboDec = comm.TurboDecoder(... 
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'NumIterations', 50, ...
    'Algorithm', 'Max*', ...
    'NumScalingBits', 3, ...
    'InputIndicesSource', 'Property', ...
    'InputIndices', outIndices);

%% 定义 BPSK 调制映射函数
% 映射规则: 0 映射为 +1, 1 映射为 -1
bpskMod = @(bits) 1 - 2*double(bits);

%% 初始化 BER 储存向量
ber_ldpc  = zeros(size(snr_dB));
ber_turbo = zeros(size(snr_dB));

%% 仿真循环
fprintf('开始仿真（突发噪声环境）：LDPC 与 Turbo 码对比\n');
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    numErrors_ldpc  = 0;  numTotal_ldpc = 0;
    numErrors_turbo = 0;  numTotal_turbo = 0;
    
    % 对于 BPSK，每维噪声方差 noiseVar = 1/(2*10^(snr/10))
    noiseVar = 1/(2*10^(snr/10));
    
    for frame = 1:numFrames
        %% 生成随机信息比特 (长度为 324)
        infoBits = randi([0 1], k, 1) > 0;
        
        %% LDPC 编码、调制及突发噪声信道传输
        % LDPC 编码
        codeword_ldpc = ldpcEncode(infoBits, encoderConfig_ldpc);
        % BPSK 调制
        txSymbols_ldpc = bpskMod(codeword_ldpc);
        
        % 生成整个码字的 AWGN 噪声
        noise_ldpc = sqrt(noiseVar) * randn(N_ldpc,1);
        % 随机判断是否发生突发噪声
        if rand < burstProb
            burstLength = round(burstFrac * N_ldpc); % 突发噪声长度
            startIdx = randi([1, N_ldpc - burstLength + 1]); % 随机突发起始位置
            % 在突发区间内，噪声方差放大 burstFactor 倍
            noise_ldpc(startIdx:startIdx+burstLength-1) = sqrt(burstFactor * noiseVar) * randn(burstLength,1);
        end
        
        % 接收信号 = 发射信号 + 噪声
        rxSymbols_ldpc = txSymbols_ldpc + noise_ldpc;
        % 计算 LDPC 部分 LLR，公式 2*y/σ² (这里 σ² 为 noiseVar)
        rxLLR_ldpc = 2 * rxSymbols_ldpc ./ noiseVar;
        % LDPC 译码
        decodedBits_ldpc = ldpcDecode(rxLLR_ldpc, decoderConfig_ldpc, maxLDPCIter, ...
                                      'OutputFormat', 'info', 'DecisionType', 'hard');
        numErrors_ldpc = numErrors_ldpc + sum(infoBits ~= decodedBits_ldpc);
        numTotal_ldpc = numTotal_ldpc + k;
        
        %% Turbo 编码、调制及突发噪声信道传输
        % Turbo 编码 (输出长度应为 648)
        encoded_turbo = step(turboEnc, infoBits);
        % BPSK 调制
        txSymbols_turbo = bpskMod(encoded_turbo);
        
        % 生成整个码字的 AWGN 噪声 (与 turbo 输出长度相同)
        noise_turbo = sqrt(noiseVar) * randn(length(encoded_turbo),1);
        if rand < burstProb
            burstLength = round(burstFrac * length(encoded_turbo));
            startIdx = randi([1, length(encoded_turbo) - burstLength + 1]);
            noise_turbo(startIdx:startIdx+burstLength-1) = sqrt(burstFactor * noiseVar) * randn(burstLength,1);
        end
        
        rxSymbols_turbo = txSymbols_turbo + noise_turbo;
        % 计算 Turbo 译码所需的 LLR 值，公式为 -2*y/σ² (注意符号及比例因子)
        rxLLR_turbo = (-2/(noiseVar)) * rxSymbols_turbo;
        % Turbo 译码
        decodedBits_turbo = step(turboDec, rxLLR_turbo);
        numErrors_turbo = numErrors_turbo + sum(infoBits ~= decodedBits_turbo);
        numTotal_turbo = numTotal_turbo + k;
    end
    ber_ldpc(s)  = numErrors_ldpc / numTotal_ldpc;
    ber_turbo(s) = numErrors_turbo / numTotal_turbo;
    
    fprintf('SNR = %.1f dB: LDPC BER = %e, Turbo BER = %e\n', snr, ber_ldpc(s), ber_turbo(s));
end

%% 绘制 BER-SNR 性能对比图
figure;
semilogy(snr_dB, ber_ldpc, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(snr_dB, ber_turbo, 'r-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('比特错误率 (BER)');
title('LDPC vs Turbo (n = 648, rate = 0.5) 在突发噪声环境下的性能');
legend('LDPC', 'Turbo');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 辅助函数：生成 H1 矩阵
function H1 = generateH1(numRows, numCols, d_v)
% generateH1: 生成一个大小为 (numRows x numCols) 的稀疏逻辑型矩阵，
% 每列恰好含有 d_v 个 1，其他元素为 0。
H1 = false(numRows, numCols);
for j = 1:numCols
    rows = randperm(numRows, d_v);
    H1(rows, j) = true;
end
H1 = sparse(H1);
end