%% ================================================================
%  1초 저항(dR) 계산 + 파일별 저장 + 전체 요약
%  - MeanBySOC: R_1s_mean, n, R_1s_std 추가
% ================================================================
clc; clear; close all;

% --- 입력/저장 경로 ---
folder_SIM = 'G:\공유 드라이브\BSL_Data4\HNE_SOC_moving_cutoff_5_processed\SIM_parsed\이름 정렬';
save_path  = 'G:\공유 드라이브\BSL_Data4\HNE_SOC_moving_cutoff_5_processed\SIM_parsed\1초 저항';
if ~exist(save_path,'dir'), mkdir(save_path); end

plot_outdir = fullfile(save_path, 'scatter_plots');
if ~exist(plot_outdir,'dir'), mkdir(plot_outdir); end

files = dir(fullfile(folder_SIM,'*_SIM.mat'));
assert(~isempty(files), 'SIM 파일을 찾지 못했습니다.');

% --- 설정 ---
soc_targets = [90 70 50 30];
soc_tol     = 1.0;   % ±1 %

% --- 전체 요약 누적 컨테이너 ---
G_file = strings(0,1);
G_sim  = zeros(0,1);
G_soc  = zeros(0,1);
G_R1s  = zeros(0,1);
G_R2   = zeros(0,1);

for f = 1:numel(files)
    F = fullfile(files(f).folder, files(f).name);
    S = load(F,'SIM_table');
    if ~isfield(S,'SIM_table'), warning('SIM_table 없음: %s', files(f).name); continue; end
    T = S.SIM_table;

    base_raw = erase(files(f).name,'_SIM.mat');

    % === 파일별 누적 컨테이너 ===
    SIM_idx   = [];
    SOC_tgt   = [];
    dt_cell   = {}; dI_cell = {}; dV_cell = {};
    OCV1_all  = []; OCV2_all = [];
    SOC1_all  = []; SOC2_all = [];
    dR_all    = []; R2_all   = [];

    for s = 1:height(T)
        SOC1 = T.SOC1(s);  SOC2 = T.SOC2(s);
        hit  = abs(SOC1 - soc_targets) <= soc_tol & abs(SOC2 - soc_targets) <= soc_tol;
        if ~any(hit), continue; end
        thisSOC = soc_targets(hit);               % 매칭된 대표 SOC(단 하나)

        % 신호
        t = T.time{s}; I = T.current{s}; V = T.voltage{s};
        if isduration(t), t = seconds(t); end
        t = t(:); I = I(:); V = V(:);
        if numel(t) < 2 || any(~isfinite(t)|~isfinite(I)|~isfinite(V)), continue; end

        % === 정수초 그리드 재샘플 → 1초 차분 ===
        t_rel = t - t(1);
        sec_marks = 0:floor(t_rel(end));
        idx_sec = arrayfun(@(ss) find(t_rel >= ss,1,'first'), sec_marks);
        idx_sec = unique(idx_sec,'stable');

        t1 = t(idx_sec);  I1 = I(idx_sec);  V1 = V(idx_sec);
        dI = diff(I1);    dV = diff(V1);    dt = diff(t1);

        ok = isfinite(dI)&isfinite(dV)&isfinite(dt)&dt>0&dt<=1.5;
        dI = dI(ok); dV = dV(ok); dt = dt(ok);

        % === 1초 저항: y = a x (절편 0) ===
        x = dI; y = dV;
        denom = sum(x.^2);
        if denom <= eps
            a = NaN; R2 = NaN;
        else
            a = sum(x.*y) / denom;                    % R_1s
            yhat = a*x;
            ss_res = sum((y - yhat).^2);
            ss_tot = sum((y - mean(y)).^2);
            R2 = 1 - ss_res/max(ss_tot, eps);
        end

        % OCV 복사(없으면 NaN)
        O1 = NaN; O2 = NaN;
        if ismember('OCV1',T.Properties.VariableNames), O1 = T.OCV1(s); end
        if ismember('OCV2',T.Properties.VariableNames), O2 = T.OCV2(s); end

        % 파일별 누적
        SIM_idx(end+1,1) = s;
        SOC_tgt(end+1,1) = thisSOC;
        dt_cell{end+1,1} = dt(:);
        dI_cell{end+1,1} = x(:);
        dV_cell{end+1,1} = y(:);
        OCV1_all(end+1,1)= O1;
        OCV2_all(end+1,1)= O2;
        SOC1_all(end+1,1)= SOC1;
        SOC2_all(end+1,1)= SOC2;
        dR_all(end+1,1)  = a;
        R2_all(end+1,1)  = R2;

        % 전체 요약 누적
        G_file(end+1,1)  = string(base_raw);
        G_sim (end+1,1)  = s;
        G_soc (end+1,1)  = thisSOC;
        G_R1s(end+1,1)   = a;
        G_R2 (end+1,1)   = R2;

        % ---- 산점도 + y=ax(절편0) 추세선 저장 ----
        if ~isempty(x) && ~isempty(y) && isfinite(a)
            % 원점-기준 R^2 (엑셀 intercept=0와 해석 일치)
            yhat   = a*x;
            SSE    = sum((y - yhat).^2);
            R2zero = 1 - SSE / max(sum(y.^2), eps);
            R2mean = R2;  % 기존 계산값 유지(평균 기준 R^2)

            fig = figure('Visible','off','Position',[100 100 720 480]);
            scatter(x, y, 18, 'filled'); grid on; hold on;
            xlabel('\Delta I (A)'); ylabel('\Delta V (V)');
            title(sprintf('%s | SIM %d | SOC %g', base_raw, s, thisSOC), 'Interpreter','none');

            % --- 축에 반드시 0 포함(시각적으로 y절편=0 보이도록) ---
            xpad = 0.05*max(1, range(x));
            ypad = 0.05*max(1, range(y));
            xlim([min([x;0]) - xpad, max([x;0]) + xpad]);
            ylim([min([y;0]) - ypad, max([y;0]) + ypad]);

            % 원점 보조선(선택)
            xline(0,'k:');
            yline(0,'k:');

            % --- y = a x 추세선: 현재 xlim 전 구간(→ 원점 포함 보장) ---
            xl = xlim;
            xx = linspace(xl(1), xl(2), 200);
            plot(xx, a*xx, ':', 'LineWidth', 2);

            % 주석(표시할 R^2 선택: R2zero 또는 R2mean)
            ax = gca;
            xpos = ax.XLim(1) + 0.05*range(ax.XLim);
            ypos = ax.YLim(2) - 0.10*range(ax.YLim);
            txt = sprintf('y = %.4g x\nR^2 = %.4f', a, R2zero);
            text(xpos, ypos, txt, 'BackgroundColor','w', 'Margin',4);

            % 파일명: <base>_SIM##_SOC##.png
            png_name = sprintf('%s_SIM%02d_SOC%02d.png', base_raw, s, round(thisSOC));
            saveas(fig, fullfile(plot_outdir, png_name));
            close(fig);
        end

    end

    % === 파일별 Summary 테이블 ===
    Summary = table( ...
        SIM_idx, SOC_tgt, ...
        dt_cell, dI_cell, dV_cell, ...
        OCV1_all, OCV2_all, SOC1_all, SOC2_all, ...
        dR_all, R2_all, ...
        'VariableNames', {'SIM','SOC_target','dt','dI','dV','OCV1','OCV2','SOC1','SOC2','R_1s','R2'});

    % === 파일별 대표 SOC 평균 + n + std ===
    soc_list  = soc_targets(:);
    mean_vals = nan(numel(soc_list),1);
    n_vals    = zeros(numel(soc_list),1);
    std_vals  = nan(numel(soc_list),1);
    for k = 1:numel(soc_list)
        idx = (SOC_tgt == soc_list(k)) & isfinite(dR_all);
        n_vals(k)    = sum(idx);
        if n_vals(k) > 0
            mean_vals(k) = mean(dR_all(idx), 'omitnan');
            std_vals(k)  = std (dR_all(idx), 0, 'omitnan');
        end
    end
    MeanBySOC = table(soc_list, n_vals, mean_vals, std_vals, ...
        'VariableNames', {'SOC_target','n','R_1s_mean','R_1s_std'});

    % 저장: <base>_R1s.mat  (→ save_path)
    out_file = fullfile(save_path, sprintf('%s_R1s.mat', base_raw));
    save(out_file, 'Summary','MeanBySOC','-v7.3');
    fprintf('[saved] %s (rows=%d)\n', out_file, height(Summary));
end



% === 전체 요약 저장 ===
Summary_all = table(G_file, G_sim, G_soc, G_R1s, G_R2, ...
    'VariableNames', {'file','SIM','SOC_target','R_1s','R2'});
save(fullfile(save_path, 'R1s_summary_all.mat'), 'Summary_all', '-v7.3');
disp('완료: R1s_summary_all.mat 저장됨.');
