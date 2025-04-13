%% simulate_ldpc_vs_turbo.m
% 本脚本对比 LDPC 和传统 Turbo 码的性能，两种码制均为：
%   码长 n = 648, 信息位数 k = 324, 码率 = 0.5
%
% LDPC 部分：使用系统形式的奇偶校验矩阵 H = [H1 | I]；
% Turbo 部分：采用 MATLAB 自带 comm.TurboEncoder 和 comm.TurboDecoder，
%   通过自定义 OutputIndices (对于编码器) 和 InputIndices (对于译码器)
%   实现 puncturing 达到率 1/2。
%
% 仿真采用 BPSK 调制、AWGN 信道，并在一定 SNR 范围内统计 BER。

%% 仿真参数
snr_dB    = 0:1:10;     % SNR (dB) 范围
numFrames = 200;        % 每个 SNR 下仿真帧数
k         = 324;        % 信息位数
N_ldpc    = 648;        % LDPC 码字长度  (rate = 324/648 = 0.5)

%% LDPC 编码器构造
M_ldpc = N_ldpc - k;    % 奇偶校验矩阵行数
d_v = 3;                % 每列1的个数
H1 = generateH1(M_ldpc, M_ldpc, d_v);       % 生成 M x M 随机矩阵 (稀疏逻辑型)
I_part = sparse(logical(speye(M_ldpc)));     % 单位矩阵 (稀疏逻辑型)
H_ldpc = [H1, I_part];      % 系统型奇偶校验矩阵

% 创建 LDPC 编码器和译码器配置对象
encoderConfig_ldpc = ldpcEncoderConfig(H_ldpc);
decoderConfig_ldpc = ldpcDecoderConfig(encoderConfig_ldpc);
maxLDPCIter = 50;  % 固定译码迭代次数

%% Turbo 码构造（采用内置 TurboEncoder/Decoder）
% 使用默认 Trellis: poly2trellis(4, [13 15], 13)
trellisTurbo = poly2trellis(4, [13 15], 13);
pMLen = log2(trellisTurbo.numStates);  % 应为 3
pN = log2(trellisTurbo.numOutputSymbols);  % 应为 2
fullOutputLen = (k + pMLen) * pN * 2;  % 全输出长度，例如 (324+3)*2*2 = 1308

% 为实现率 1/2 (即传输 648 个符号)，设计简单均匀的 puncturing 索引：
% 从 fullOutputLen=1308 个符号中，取前 1296 个，然后均匀选择 648 个符号
tempIdx = (1:2:1296).';    % 产生 648 个索引
outIndices = tempIdx;      % 用作 TurboEncoder 的 OutputIndices
% 对于 TurboDecoder，使用相同的 puncturing pattern，赋值给 InputIndices

% 设定简单的交织器，这里采用逆序作为例子
interleaverIndices = (k:-1:1).';

% 创建 Turbo 编码器
turboEnc = comm.TurboEncoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'OutputIndicesSource', 'Property', ...
    'OutputIndices', outIndices);

% 创建 Turbo 译码器，注意使用 InputIndices 属性指定 puncturing 信息
turboDec = comm.TurboDecoder(...
    'TrellisStructure', trellisTurbo, ...
    'InterleaverIndices', interleaverIndices, ...
    'NumIterations', 6, ...           % 可根据需要调整译码迭代次数
    'InputIndicesSource', 'Property', ...
    'InputIndices', outIndices);

%% 定义 BPSK 调制映射函数
bpskMod = @(bits) 1 - 2*double(bits);

%% 初始化 BER 储存矩阵
ber_ldpc   = zeros(size(snr_dB));
ber_turbo  = zeros(size(snr_dB));

%% 仿真
fprintf('开始仿真：LDPC 与 Turbo 码对比\n');
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    numErrors_ldpc = 0;  numTotal_ldpc = 0;
    numErrors_turbo = 0; numTotal_turbo = 0;
    
    % 对于 BPSK, 噪声方差 noiseVar = 1/(2*10^(snr/10))
    noiseVar = 1/(2*10^(snr/10));
    
    for frame = 1:numFrames
        %% 生成随机信息比特 (324x1)
        infoBits = randi([0 1], k, 1) > 0;
        
        %% LDPC 编码与译码
        codeword_ldpc = ldpcEncode(infoBits, encoderConfig_ldpc);
        txSymbols_ldpc = bpskMod(codeword_ldpc);
        rxSymbols_ldpc = txSymbols_ldpc + sqrt(noiseVar)*randn(N_ldpc,1);
        rxLLR_ldpc = 2*rxSymbols_ldpc./noiseVar;
        decodedBits_ldpc = ldpcDecode(rxLLR_ldpc, decoderConfig_ldpc, maxLDPCIter, ...
                          'OutputFormat', 'info', 'DecisionType', 'hard');
        numErrors_ldpc = numErrors_ldpc + sum(infoBits ~= decodedBits_ldpc);
        numTotal_ldpc = numTotal_ldpc + k;
        
        %% Turbo 编码与译码
        encoded_turbo = step(turboEnc, infoBits);  % 输出长度应为 648
        txSymbols_turbo = bpskMod(encoded_turbo);
        rxSymbols_turbo = txSymbols_turbo + sqrt(noiseVar)*randn(length(encoded_turbo),1);
        % 对于 Turbo 解码器，按照示例计算 LLR
        rxLLR_turbo = (-2/(noiseVar/2)) * real(rxSymbols_turbo);
        decodedBits_turbo = step(turboDec, rxLLR_turbo);
        numErrors_turbo = numErrors_turbo + sum(infoBits ~= decodedBits_turbo);
        numTotal_turbo = numTotal_turbo + k;
    end
    ber_ldpc(s)   = numErrors_ldpc / numTotal_ldpc;
    ber_turbo(s)  = numErrors_turbo / numTotal_turbo;
    
    fprintf('SNR = %.1f dB: LDPC BER = %e, Turbo BER = %e\n', snr, ber_ldpc(s), ber_turbo(s));
end

%% 绘制对比图
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
%% 辅助函数：生成 H1 部分 (numRows x numCols 稀疏逻辑型矩阵，每列 d_v 个1)
function H1 = generateH1(numRows, numCols, d_v)
% generateH1 生成一个 numRows x numCols 的稀疏逻辑型矩阵，
% 每列包含 d_v 个1，其余为0.
H1 = false(numRows, numCols);
for j = 1:numCols
    rows = randperm(numRows, d_v);
    H1(rows, j) = true;
end
H1 = sparse(logical(H1));
end