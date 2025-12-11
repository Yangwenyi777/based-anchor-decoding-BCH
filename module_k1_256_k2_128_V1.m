%% ----------代码思路------------
% 1、首先随机生成4KB的数据
% 2、将4KB数据转换成256*128的矩阵形式
% 3、计算256的BCH纠错码，要求能够纠正2个错误(由于BCH对码长有要求，因此需要缩短码)
% 4、计算128的BCH码，要求能够纠正3个错误
% 5、将新计算好的272*134矩阵送入BSC信道进行随机注错
% 6、对注错后的矩阵进行解码，之后与注错前的矩阵比较，记录错误bit和纠错能力
% 7、回到第1点重复上述操作，直到超过最大错误bit数或者纠错能力达到目标值跳出循环
% 8、画出对应的曲线图
% 
clc; clear;close all;
% if isempty(gcp('nocreate'))
%     parpool('local', 8);   % i7-11800：8 个 worker
% end
%% ---------- 参数配置 ----------
dataBytes      = 4096;                  % 4KB 数据
bitsPerByte    = 8;
% targetErrBits  = 1e3;                   % 终止条件之一：累计错比特
% maxBits        = 1e8;                   % 终止条件之二：最大发送比特数量
useAWGN        = false;                 % 选择 AWGN 注错; false = 手动注错
max_loop       = 10;                    % 迭代纠错的循环次数
max_test_loop  = 27300;
% 1e7的数据量是max_test_loop = 273;  mins
% 1e8的数据量是max_test_loop = 2730; hours
RBER = 0.015 : -0.001 : 0.005;            % 在 10^-1=0.1 到 10^-3 之间生成13个点
% RBER = logspace(-1, -3, 10);      
% 几何：4KB -> 256x128 比特矩阵
nRows_info     = 256;                 % 2)
nCols_info     = 128;

% 行向 BCH：采用合法的 BCH(274,256,t=2) 利用RS(511,493)进行编码，不足的k用0补全
n_bch_full     = 274; 
k_bch_full     = 256;
n_row_bits     = 511;
k_row_bits     = 493;
shorten_bch    = k_row_bits - k_bch_full;   % 493 - 256= 237

% 列向 BCH：按"字节列"做 RS(152,128,t=3)（每符号=1字节）
N_rs            = 144; 
K_rs            = 128;                   % t=(152-128)/8 = 3
m_rs            = 8;                     % 1 个符号 = 8 bit
n_col_bits      = 255;
k_col_bits      = 239; 
shorten_col_bch = k_col_bits - K_rs;

% 码率：行率 * 列率
rate_row       = k_bch_full / n_bch_full;   % 128/144
rate_col       = K_rs / N_rs;               % 128/134
rate_sys       = rate_row * rate_col;

% 误码统计
errBER         = zeros(numel(RBER),1);   %用一个数组存储变量BER(bit error rate)
correct_rate   = zeros(numel(RBER), max_test_loop);

debug_printf   = false;
%% ---------- 主循环 ----------
for ie = 1:numel(RBER)
    
    hBER = comm.ErrorRate;   
    col_corate = int32(1);
    error_num = 0;
    error_num_before = 0;
    for test_num = 1:max_test_loop
%         fprintf('开始时间：%s\n', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS')));
        % 1) 随机生成 4KB 数据
        txBitsAll = randi([0 1], dataBytes*bitsPerByte, 1);
%         rng(1234);   % 固定随机种子
%         txBitsAll = uint8(randi([0 1], 32768, 1));

        % 2) 4KB -> 256x128 比特矩阵（逐行填充）
        infoMat = reshape(txBitsAll, nCols_info, nRows_info);  % 256x128        

        % 3) 列向 BCH 编码：按"字节列"做 BCH(144,128,t=3)
%         extend_col_bch_matrix = false(n_col_bits, k_bch_full);
        temp1 = infoMat.';
        temp2 = zeros(k_bch_full, shorten_col_bch, 'double');
        msg240 = [temp2, temp1];
        gf_msg240 = gf(msg240, 1);
        cw256 = bchenc(gf_msg240, n_col_bits, k_col_bits);
        double_cw256 = double(cw256.x);
        extend_col_bch_matrix = double_cw256;
        col_end_matrix = extend_col_bch_matrix(:,shorten_col_bch+1:n_col_bits);

        % 4) 行向 BCH 编码：BCH(511,493)
        %    做法：每行前补 237 个 0 → 493位，BCH 编码得 511 位
%         rowCwBits = false(n_row_bits, N_rs).';
        msg493    = [zeros(N_rs, shorten_bch, 'double'), col_end_matrix.'];
        gf_msg493 = gf(msg493, 1);
        cw511     = bchenc(gf_msg493, n_row_bits, k_row_bits);
        tx_cw511 = double(cw511.x);
        rowCwBits = tx_cw511;
        eBCH2_matrix = logical(rowCwBits(:,shorten_bch+1:n_row_bits));

%         temp4 = eBCH2_matrix(1:128,257:274).';
%         temp_msg = [zeros(18, shorten_col_bch, 'double'), temp4];
%         temp_gf_msg = gf(temp_msg, 1);
%         temp_gf_msg1 = bchenc(temp_gf_msg, 255, 239);
%         temp_msg1 = double(temp_gf_msg1.x);
%         temp5 = temp_msg1(: , shorten_col_bch+1:n_col_bits);
%         temp6 = eBCH2_matrix(1:144,257:274);
%         diff_mat = (temp5.' ~= temp6);
%         [r,c] = find(diff_mat);



        % 5) 通过信道：AWGN 或手动注错
        txBits = eBCH2_matrix(:);
        if useAWGN
            txSig  = qammod(txBits, M, 'UnitAveragePower', true, 'InputType','bit');
            rxSig  = awgn(txSig, SNRdB);
            rxBits = qamdemod(rxSig, M, 'UnitAveragePower', true, 'OutputType','bit');
        else
            rxBits = txBits;
            % 手动随机翻转比特（你可以自定义翻转位置/数量）
            flipNum = max(1, round(RBER(1,ie)*numel(rxBits)));
            flipIdx = randperm(numel(rxBits), flipNum);
            rxBits(flipIdx) = ~rxBits(flipIdx);
        end

        % 6) 译码:先译码行,后译码列
        rx_infoMat_144x274 = reshape(rxBits, N_rs, n_bch_full);          % 134x274
        if debug_printf
            fprintf('比较tx和rx信号,对比注错过程\n');
            fprintf('原始误码率=%d',RBER(1,ie));
            report_info_changes(eBCH2_matrix,rx_infoMat_144x274,N_rs,274);   % = flipNum, judge error num in agwn, should be same with flipNum
        end
        for j = 1:max_loop
            % 逆 BCH：基于274bit的基础上补全511bit的码长，送入解码器
            backup_rxMat = rx_infoMat_144x274;   
            % 一次性构造 144×511 的矩阵
            rxRowSrc = [zeros(N_rs, shorten_bch, 'logical'), rx_infoMat_144x274];            
            % 转成 GF(2)，MATLAB 要求每行是一个码字
            gf_rxRowSrc = gf(rxRowSrc, 1);            
            % 一次性对 144 个码字进行 BCH 解码
            [~, ~, gf_rx_msg511] = bchdec(gf_rxRowSrc, n_row_bits, k_row_bits);           
            % 取回 274 bit 信息部分（对所有行一次性完成）
            rx_msg511 = double(gf_rx_msg511.x);
            rx_infoMat_144x274 = rx_msg511(:, shorten_bch+1 : n_row_bits);  
            if debug_printf
                fprintf('行BCH译码完成后,比较一下纠正了多少错误');
                report_info_changes(backup_rxMat, rx_infoMat_144x274, N_rs, 274);
                fprintf('行BCH译码完成后,比较一下还剩多少错误');
                report_info_changes(eBCH2_matrix,rx_infoMat_144x274,N_rs,274);
            end
     
            % -------- 逆 BCH（列方向）-----------
            beforeRS_mat = rx_infoMat_144x274;
            % rx_infoMat_144x274: 144×274
            % 先把它转成 274×144，再左侧补零 shorten_col_bch 到 n_col_bits=255
            rxColSrc = [zeros(shorten_col_bch, n_bch_full,'logical').', rx_infoMat_144x274.'];
            % 转 GF
            gf_rxColSrc = gf(rxColSrc, 1);  
            % 一次性对 274 个码字进行 BCH 解码
            [~, ~, gf_decoder_msg255] = bchdec(gf_rxColSrc, n_col_bits, k_col_bits);  
            % 转为 double，取后 144 bit 信息区
            decoder_msg255 = double(gf_decoder_msg255.x);  
            % 写回（注意矩阵维度）
            rx_infoMat_144x274 = decoder_msg255(:, shorten_col_bch+1 : n_col_bits).';
            if debug_printf
                fprintf('列BCH译码完成后,比较一下列纠正了多少错误');
                report_info_changes(beforeRS_mat,rx_infoMat_144x274,N_rs,256);       %debug
                fprintf('列BCH译码完成后,比较一下还剩多少错误');
                report_info_changes(eBCH2_matrix,rx_infoMat_144x274,N_rs,274);
            end
            temp3 = logical(rx_infoMat_144x274);
            [errmatrix, num_changed_total] = report_info_changes(backup_rxMat,temp3,N_rs,274);       %debug
            if(num_changed_total == 0)
                if debug_printf
                    fprintf('3333333:');
                    report_info_changes(eBCH2_matrix,rx_infoMat_144x274,N_rs,274);
                end
                break;
            end
        end
        
        txInfoBits = eBCH2_matrix(:);
        rxInfoBits = logical(rx_infoMat_144x274(:));

        stats = hBER(txInfoBits, rxInfoBits);
        error_num = stats(2) - error_num_before;
        correct_num = flipNum - error_num;
        correct_rate(ie, max_test_loop) = correct_num / flipNum;
%         col_corate = col_corate + 1;
        error_num_before = stats(2);
        fprintf('cur_Errs=%d correct_errs = %d\n',stats(2), correct_num);
%         if stats(2) >= targetErrBits || stats(3) >= maxBits
%             break;
%         end
%         fprintf('结束时间：%s\n\n', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS')));
    end

    errBER(ie) = stats(1);      %BER = stats(1)
    fprintf('RPBER=%.1f dB | BER=%.3e | send_Bits=%d | Errs=%d\n', ...
        RBER(ie), errBER(ie), stats(3), stats(2));
    
end

%% ---------- 作图 ----------
semilogy(RBER, errBER, 'o-'); grid on;
xlabel('RBER (dB)'); ylabel('BER');
title('TPC: row BCH(t=2, 256->274) + column RS(134,128,t=3)');
legend('Sim BER','Location','southwest');
