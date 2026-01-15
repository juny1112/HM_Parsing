clc; clear; close all;

% ===== 처리할 MAT 파일 "하나" 경로 =====
mat_path = "G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\parsed_data\HNE_10degC_1C_13_7_1113.mat";

% ===== Fig 저장 폴더(원하는 곳) =====
fig_folder = "C:\Users\junny\Downloads";
if ~exist(fig_folder, 'dir'), mkdir(fig_folder); end

S = load(mat_path);

% ---- 데이터 꺼내기 -------------------------------------------------
if isfield(S, 'parsed_data')
    pd = S.parsed_data;

    valid = arrayfun(@(s) ~isempty(s.time), pd);
    pd = pd(valid);
    if isempty(pd)
        warning('유효한 스텝이 없습니다. 건너뜀: %s', mat_path);
        return;
    end

    t_all = vertcat(pd.time);      % duration
    v_all = vertcat(pd.voltage);
    i_all = vertcat(pd.current);

elseif isfield(S, 'intrim')
    t_all = S.intrim.time;         % duration
    v_all = S.intrim.voltage;
    i_all = S.intrim.current;

else
    warning('parsed_data/intrim 모두 없음. 건너뜀: %s', mat_path);
    return;
end

t_all = t_all(:); v_all = v_all(:); i_all = i_all(:);
t_sec = seconds(t_all - t_all(1));

% ---- 플롯 ----------------------------------------------------------
fig = figure('Visible','on');
yyaxis left
plot(t_sec, i_all, 'b-', 'LineWidth', 1); hold on
ylabel('Current [A]')

yyaxis right
plot(t_sec, v_all, 'r-', 'LineWidth', 1); hold off
ylabel('Voltage [V]')

[~, base, ~] = fileparts(mat_path);
title(strrep(base, '_', '\_'), 'Interpreter','tex')
xlabel('Time [s]')
box on; grid on;
set(gcf, 'Position', [300 300 1200 800])

% ---- 저장 ----------------------------------------------------------
fig_name = fullfile(fig_folder, base + ".fig");
savefig(fig, fig_name);

fprintf("저장 완료: %s\n", fig_name);
