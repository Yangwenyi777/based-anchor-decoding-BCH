function [delta_data, num_changed_data, pos, row_bit_changes, row_err_num] = report_info_changes(infoMat, rxMat_134x274, K, C)
% REPORT_INFO_CHANGES
% 对比信息区 (前 K 行、前 C 列) 的比特是否变化，并打印变化的比特总数。
% 默认 K=128, C=256（与你当前代码一致）。
%
% 返回：
%   delta_data        - 逻辑矩阵，1 表示该比特与原始不同
%   num_changed_data  - 变化的比特总数

    if nargin < 3, K = 128; end
    if nargin < 4, C = 256; end

    % 取出信息区，并转换为逻辑类型
    A = logical(infoMat(1:K, 1:C));
    B = logical(rxMat_134x274(1:K, 1:C));

    % bit 级变化
    delta_data = xor(A, B);
    num_changed_data = nnz(delta_data);
    row_bit_changes = sum(delta_data, 2);   
    row_err_num = nnz(row_bit_changes);
    [r, c] = find(delta_data);
    pos = [r, c];
    % 打印
    fprintf('[CHG[data]] changed info bits = %d, changed symbol number = %d  (region %dx%d)\n', ...
        num_changed_data, row_err_num, K, C);
end