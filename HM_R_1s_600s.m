%% ================================================================
%  대표 SOC(±tol) 풀링 기반 1초 저항(dR) 계산 — SOC_vec 기준, 보간 없음
%  - 정수초 그리드에서 '첫 샘플' 선택 → 1초 차분
%  - (옵션) SOC_vec ∈ [SOC±tol]인 페어만 사용 (tol 미설정이면 SOC 필터 없이 전부 사용)
%  - 원점회귀
%  - (NEW) 각 SIM 세그는 0~600초 구간만 사용
%  - (NEW) US06 뿐만 아니라 모든 주행부하(8종) 각각에 대해 정리
%  - (ADD) do_plot 스위치로 (file,load,SOC)별 pooled scatter + 회귀선 PNG 저장
% ================================================================
clc; clear; close all;

%% 경로/설정 -------------------------------------------------------------
folder_SIM = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\10degC\이름정렬';

fit_window_sec = 600;  % ★ 600초 제한

save_path  = fullfile(folder_SIM, sprintf('1초저항_allLoads_%ds', fit_window_sec));

% ── 산점도 plot on/off --------------------------------------------------
do_plot = true;   % ← 산점도 저장하고 싶으면 true

if ~exist(save_path,'dir'), mkdir(save_path); end
plot_outdir = fullfile(save_path, 'scatter_pooled');
if do_plot && ~exist(plot_outdir,'dir'), mkdir(plot_outdir); end

files = dir(fullfile(folder_SIM,'*_SIM.mat'));
assert(~isempty(files), 'SIM 파일을 찾지 못했습니다.');

% 대표 SOC
soc_targets = [50 70];

% ★ SOC 허용오차(옵션)
%   - [] 또는 NaN: SOC 필터 없이 0~600초 전체 페어 사용
%   - 숫자: abs(SOC_mid - tgt) <= soc_tol 적용
soc_tol = [];   % 예) 5 로 주면 ±5% 필터, []면 필터 없음

% 주행부하 순서
loadNames = {'US06','UDDS','HWFET','WLTP','CITY1','CITY2','HW1','HW2'};
nLoads    = numel(loadNames);

% 블록 크기 = 부하 개수
blkSize = nLoads;

%% 전체 요약 누적 (R1은 mΩ)
G_file = strings(0,1);
G_load = strings(0,1);
G_soc  = [];
G_R1s  = [];
G_R2   = [];
G_npts = [];
G_SIM  = strings(0,1);

%% 파일 루프 -------------------------------------------------------------
for f = 1:numel(files)
    F = fullfile(files(f).folder, files(f).name);
    base_raw = erase(files(f).name,'_SIM.mat');

    S = load(F,'SIM_table');
    if ~isfield(S,'SIM_table')
        warning('SIM_table 없음: %s', files(f).name);
        continue;
    end
    T = S.SIM_table;

    % 필수 컬럼 확인
    needCols = {'time','current','voltage','SOC_vec'};
    if ~all(ismember(needCols, T.Properties.VariableNames))
        warning('[%s] 필수 컬럼 누락: %s', base_raw, ...
            strjoin(setdiff(needCols, T.Properties.VariableNames), ', '));
        continue;
    end

    hasRowNames = ~isempty(T.Properties.RowNames);

    K = numel(soc_targets);
    Result_byLoad = struct;

    for l = 1:nLoads
        loadName = loadNames{l};

        pooled = repmat(struct('x',[],'y',[]), K, 1);
        n_pairs = zeros(K,1);
        R1s_mOhm = nan(K,1);
        R2       = nan(K,1);
        SIM_used = strings(K,1);

        for k = 1:K
            sim_idx_expect = (k-1)*blkSize + l;

            % row 선택
            sRow = NaN;
            if hasRowNames
                sim_name = sprintf('SIM%d', sim_idx_expect);
                hit = find(strcmp(T.Properties.RowNames, sim_name), 1, 'first');
                if ~isempty(hit)
                    sRow = hit;
                    SIM_used(k) = string(sim_name);
                end
            end

            % fallback: RowNames 없거나 못 찾으면 row index 사용
            if isnan(sRow)
                if sim_idx_expect >= 1 && sim_idx_expect <= height(T)
                    sRow = sim_idx_expect;
                    SIM_used(k) = sprintf('ROW%d', sim_idx_expect);
                else
                    SIM_used(k) = "NaN";
                    continue;
                end
            end

            % 데이터 꺼내기
            t  = T.time{sRow};
            I  = T.current{sRow};
            V  = T.voltage{sRow};
            SV = T.SOC_vec{sRow};

            if isempty(t) || numel(t) < 2, continue; end

            % 시간 → 초(double), 0 기준
            if isdatetime(t)
                t = seconds(t - t(1));
            elseif isduration(t)
                t = seconds(t); t = t - t(1);
            else
                t = t - t(1);
            end

            t = t(:); I = I(:); V = V(:); SV = SV(:);

            % 중복 time 제거
            [t, ia] = unique(t,'stable');
            I = I(ia); V = V(ia); SV = SV(ia);

            if numel(t) < 2 || any(~isfinite(t)|~isfinite(I)|~isfinite(V)|~isfinite(SV))
                continue;
            end

            % SOC 스케일 보정(0~1 → %)
            medSOC = median(SV(isfinite(SV)));
            if isfinite(medSOC) && medSOC <= 1.5, SV = SV*100; end

            % ★ 600초 트리밍
            t_rel = t - t(1);
            mWin  = (t_rel >= 0) & (t_rel <= fit_window_sec);
            if nnz(mWin) < 3
                continue;
            end
            t_rel = t_rel(mWin);
            I = I(mWin); V = V(mWin); SV = SV(mWin);

            % 정수초 그리드 재샘플(보간 아님): 각 정수초에서 "첫 샘플" 선택
            maxSec = floor(t_rel(end));
            if maxSec < 1
                continue;
            end
            sec_marks = 0:maxSec;
            idx_sec = arrayfun(@(ss) find(t_rel >= ss,1,'first'), sec_marks);
            idx_sec = unique(idx_sec,'stable');

            t1 = t_rel(idx_sec);
            I1 = I(idx_sec);
            V1 = V(idx_sec);
            S1 = SV(idx_sec);

            if numel(t1) < 2
                continue;
            end

            % 1초 차분
            dI = diff(I1);
            dV = diff(V1);
            dt = diff(t1);
            SOC_mid = 0.5*(S1(1:end-1) + S1(2:end));

            ok_base = isfinite(dI) & isfinite(dV) & isfinite(dt) & isfinite(SOC_mid) ...
                      & dt>0 & dt<=1.5;

            % ===== (핵심) SOC 필터 옵션화 =====
            if isempty(soc_tol) || (isscalar(soc_tol) && isnan(soc_tol))
                % tol 미설정: 0~600초 전체 페어 사용 (SOC 무시)
                pick = ok_base;
            else
                % tol 설정: 기존대로 SOC 근처만 사용
                tgt  = soc_targets(k);
                pick = ok_base & (abs(SOC_mid - tgt) <= soc_tol);
            end

            n_add = sum(pick);
            if n_add > 0
                pooled(k).x = [pooled(k).x; dI(pick)];
                pooled(k).y = [pooled(k).y; dV(pick)];
            end
        end

        % ---- SOC별 원점회귀 + (옵션) 산점도 저장 ----
        for k = 1:K
            x = pooled(k).x;
            y = pooled(k).y;
            n_pairs(k) = numel(x);

            if isempty(x) || sum(x.^2) <= eps
                continue;
            end

            a_ohm = sum(x.*y) / sum(x.^2);   % Ω
            yhat  = a_ohm*x;
            R2z   = 1 - sum((y - yhat).^2) / max(sum(y.^2), eps);

            R1s_mOhm(k) = a_ohm * 1000;  % mΩ
            R2(k)       = R2z;

            % ===== 산점도 저장(옵션) =====
            if do_plot
                if isempty(soc_tol) || (isscalar(soc_tol) && isnan(soc_tol))
                    tolTxt = 'tol=NONE';
                    tolTag = 'tolNONE';
                else
                    tolTxt = sprintf('tol=%.3g', soc_tol);
                    tolTag = sprintf('tol%g', soc_tol);
                end

                loadTag = regexprep(loadName, '[^\w\-]', '_');
                baseTag = regexprep(base_raw,  '[^\w\-]', '_');

                fig = figure('Visible','off','Position',[100 100 720 520]);
                scatter(x, y, 8, 'filled'); grid on; hold on;
                xlabel('\Delta I (A)'); ylabel('\Delta V (V)');

                title(sprintf('%s | %s | SOC %g pooled (n=%d) | win=%ds | %s', ...
                    base_raw, loadName, soc_targets(k), n_pairs(k), fit_window_sec, tolTxt), ...
                    'Interpreter','none');

                if numel(x) >= 2
                    xpad = 0.05*max(1, max(x)-min(x));
                    ypad = 0.05*max(1, max(y)-min(y));
                else
                    xpad = 1; ypad = 1;
                end
                xlim([min([x;0]) - xpad, max([x;0]) + xpad]);
                ylim([min([y;0]) - ypad, max([y;0]) + ypad]);
                xline(0,'k:'); yline(0,'k:');

                xl = xlim;
                xx = linspace(xl(1), xl(2), 200);
                plot(xx, a_ohm*xx, ':', 'LineWidth', 2);

                ax = gca;
                xpos = ax.XLim(1) + 0.05*range(ax.XLim);
                ypos = ax.YLim(2) - 0.18*range(ax.YLim);
                txt  = sprintf('y = %.4g x (Ω)\nR_{1s} = %.3f mΩ\nR^2 = %.4f\nn = %d', ...
                               a_ohm, R1s_mOhm(k), R2z, n_pairs(k));
                text(xpos, ypos, txt, 'BackgroundColor','w','Margin',4,'Interpreter','none');

                out_png = sprintf('%s_%s_SOC%02d_win%ds_%s.png', ...
                    baseTag, loadTag, round(soc_targets(k)), fit_window_sec, tolTag);

                saveas(fig, fullfile(plot_outdir, out_png));
                close(fig);
            end
        end

        % ---- 테이블 정리/저장(파일 내부 load별) ----
        LoadBySOC = table( ...
            repmat(string(loadName), K, 1), ...
            soc_targets(:), ...
            SIM_used(:), ...
            n_pairs(:), ...
            R1s_mOhm(:), ...
            R2(:), ...
            'VariableNames', {'Load','SOC_target','SIM_used','n_pairs','R_1s_mOhm','R2'} );

        Result_byLoad.(matlab.lang.makeValidName(loadName)) = LoadBySOC;

        % ---- 전역 요약 누적 ----
        G_file = [G_file; repmat(string(base_raw), K, 1)];
        G_load = [G_load; repmat(string(loadName), K, 1)];
        G_soc  = [G_soc;  soc_targets(:)];
        G_SIM  = [G_SIM;  SIM_used(:)];
        G_R1s  = [G_R1s;  R1s_mOhm(:)];
        G_R2   = [G_R2;   R2(:)];
        G_npts = [G_npts; n_pairs(:)];

        fprintf('[%s] %s 완료 (window=%ds)\n', base_raw, loadName, fit_window_sec);
        for k = 1:K
            if isempty(soc_tol) || (isscalar(soc_tol) && isnan(soc_tol))
                tolTxt = 'tol=NONE';
            else
                tolTxt = sprintf('tol=%.3g', soc_tol);
            end
            fprintf('   SOC %g | %s | %s | R1=%.3f mΩ | n=%d | R2=%.3f\n', ...
                soc_targets(k), SIM_used(k), tolTxt, R1s_mOhm(k), n_pairs(k), R2(k));
        end
    end

    % ---- 파일별 mat 저장 ----
    save(fullfile(save_path, sprintf('%s_R1s_pooled_allLoads_%ds.mat', base_raw, fit_window_sec)), ...
        'Result_byLoad','soc_targets','soc_tol','loadNames','fit_window_sec','-v7.3');
end

%% 전체 요약 및 피벗 저장 ----------------------------------------------
Summary_all = table(G_file, G_load, G_soc, G_SIM, G_R1s, G_R2, G_npts, ...
    'VariableNames', {'file','load','SOC_target','SIM_used','R_1s_mOhm','R2','n_pairs'});

cellNames = unique(Summary_all.file, 'stable');
socLevels = unique(Summary_all.SOC_target, 'stable');

Tbl_SOC_row_Cell_col = struct;   % load별 (행=SOC, 열=Cell)
Tbl_Cell_row_SOC_col = struct;   % load별 (행=Cell, 열=SOC)
R1_SOCxCellxLoad_mOhm = nan(numel(socLevels), numel(cellNames), nLoads);

rowNames1 = arrayfun(@(s) sprintf('SOC %d',s), socLevels, 'uni',false);
varNamesCell = matlab.lang.makeValidName(cellNames);
varNamesSOC  = arrayfun(@(s) sprintf('SOC_%d',s), socLevels, 'uni',false);

for l = 1:nLoads
    loadName = loadNames{l};
    maskL = strcmp(string(Summary_all.load), string(loadName));
    SL = Summary_all(maskL,:);

    R1_SOCxCell_mOhm = nan(numel(socLevels), numel(cellNames));

    for i = 1:height(SL)
        r = find(socLevels == SL.SOC_target(i));
        c = find(strcmp(cellNames, SL.file(i)));
        R1_SOCxCell_mOhm(r,c) = SL.R_1s_mOhm(i);
    end

    R1_SOCxCellxLoad_mOhm(:,:,l) = R1_SOCxCell_mOhm;

    Tbl1 = array2table(R1_SOCxCell_mOhm, ...
        'RowNames', rowNames1, 'VariableNames', varNamesCell);

    Tbl2 = array2table(R1_SOCxCell_mOhm.', ...
        'RowNames', varNamesCell, 'VariableNames', varNamesSOC);

    Tbl_SOC_row_Cell_col.(matlab.lang.makeValidName(loadName)) = Tbl1;
    Tbl_Cell_row_SOC_col.(matlab.lang.makeValidName(loadName)) = Tbl2;
end

save(fullfile(save_path, sprintf('R1s_pooled_allLoads_%ds_summary.mat', fit_window_sec)), ...
    'Summary_all','soc_targets','soc_tol','loadNames','fit_window_sec', ...
    'cellNames','socLevels', ...
    'R1_SOCxCellxLoad_mOhm','Tbl_SOC_row_Cell_col','Tbl_Cell_row_SOC_col','-v7.3');

disp('완료: 600초 제한 + (옵션)SOC tol 미설정 시 전체 페어 사용 + 모든 주행부하 1초저항 저장 + (옵션) 산점도 저장 완료.');
