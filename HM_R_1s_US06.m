%% ================================================================
%  대표 SOC(±tol) 풀링 기반 1초 저항(dR) 계산 — SOC_vec 기준, 보간 없음
%  - 정수초 그리드에서 '첫 샘플' 선택 → 1초 차분
%  - SOC_vec ∈ [SOC±tol]인 페어만 사용, 원점회귀
%  - SIM1, SIM9 (US06)만 사용 (RowNames 기반)
% ================================================================
clc; clear; close all;

%% 경로/설정 -------------------------------------------------------------
folder_SIM = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\20degC\이름정렬';
save_path  = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\20degC\1초저항_US06';

do_plot = false;

if ~exist(save_path,'dir'), mkdir(save_path); end
plot_outdir = fullfile(save_path, 'scatter_pooled');
if do_plot && ~exist(plot_outdir,'dir'), mkdir(plot_outdir); end

files = dir(fullfile(folder_SIM,'*_SIM.mat'));
assert(~isempty(files), 'SIM 파일을 찾지 못했습니다.');

% 대표 SOC/허용오차
soc_targets = [50 70];
soc_tol     = 5;   % ±1 %

% (중요) 사용할 SIM 이름: SIM1, SIM9
sim_pick_names = {'SIM1','SIM9'};

% 전체 요약 누적 (R1은 mΩ)
G_file = strings(0,1);
G_soc  = [];
G_R1s  = [];
G_R2   = [];
G_npts = [];

%% 파일 루프 -------------------------------------------------------------
for f = 1:numel(files)
    F = fullfile(files(f).folder, files(f).name);
    base_raw = erase(files(f).name,'_SIM.mat');

    S = load(F,'SIM_table');
    if ~isfield(S,'SIM_table'), warning('SIM_table 없음: %s', files(f).name); continue; end
    T = S.SIM_table;

    % 필수 컬럼 확인
    needCols = {'time','current','voltage','SOC_vec'};
    if ~all(ismember(needCols, T.Properties.VariableNames))
        warning('[%s] 필수 컬럼 누락: %s', base_raw, ...
            strjoin(setdiff(needCols, T.Properties.VariableNames), ', '));
        continue;
    end

    % ===== SIM1, SIM9 선택 (RowNames 기반) =====
    if ~isempty(T.Properties.RowNames)
        sim_use = find(ismember(T.Properties.RowNames, sim_pick_names));
        sim_use_names = T.Properties.RowNames(sim_use);
    else
        % RowNames가 없을 때를 대비한 방어 로직(필요시 수정)
        error('[%s] SIM_table에 RowNames가 없습니다. (SIM1/SIM9 선택 불가)', base_raw);
    end

    if isempty(sim_use)
        warning('[%s] 요청한 SIM(%s) 행이 없습니다. RowNames 확인 필요.', ...
            base_raw, strjoin(sim_pick_names, ', '));
        continue;
    end

    K = numel(soc_targets);
    pooled = repmat(struct('x',[],'y',[]), K, 1);
    counts = zeros(K, height(T));   % SOC별 × SIM별 페어 수

    % ===== SIM1, SIM9만 처리 =====
    for s = sim_use(:).'
        t  = T.time{s};
        I  = T.current{s};
        V  = T.voltage{s};
        SV = T.SOC_vec{s};
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
        [t, ia] = unique(t,'stable');
        I = I(ia); V = V(ia); SV = SV(ia);
        if numel(t) < 2 || any(~isfinite(t)|~isfinite(I)|~isfinite(V)|~isfinite(SV)), continue; end

        % SOC 스케일 보정(0~1 → %)
        medSOC = median(SV(isfinite(SV)));
        if isfinite(medSOC) && medSOC <= 1.5, SV = SV*100; end

        % 정수초 그리드 재샘플(보간 아님)
        t_rel = t - t(1);
        sec_marks = 0:floor(t_rel(end));
        idx_sec = arrayfun(@(ss) find(t_rel >= ss,1,'first'), sec_marks);
        idx_sec = unique(idx_sec,'stable');

        t1 = t(idx_sec);
        I1 = I(idx_sec);
        V1 = V(idx_sec);
        S1 = SV(idx_sec);

        % 1초 차분
        dI = diff(I1);
        dV = diff(V1);
        dt = diff(t1);
        SOC_mid = 0.5*(S1(1:end-1) + S1(2:end));

        ok_base = isfinite(dI) & isfinite(dV) & isfinite(dt) & isfinite(SOC_mid) ...
                  & dt>0 & dt<=1.5;

        % 대표 SOC별 풀링
        for k = 1:K
            tgt = soc_targets(k);
            pick = ok_base & (abs(SOC_mid - tgt) <= soc_tol);
            n_add = sum(pick);
            if n_add>0
                pooled(k).x = [pooled(k).x; dI(pick)];
                pooled(k).y = [pooled(k).y; dV(pick)];
                counts(k, s) = counts(k, s) + n_add;
            end
        end
    end

    % ===== 회귀/저장 =====
    R1s_mOhm = nan(K,1);
    R2       = nan(K,1);
    n_pairs  = zeros(K,1);
    SIM_used = cell(K,1);
    pairs_per_SIM = cell(K,1);

    for k = 1:K
        x = pooled(k).x; y = pooled(k).y;
        n_pairs(k) = numel(x);

        used_idx = find(counts(k,:) > 0);
        SIM_used{k}      = used_idx(:);
        pairs_per_SIM{k} = counts(k, used_idx).';

        if isempty(x) || sum(x.^2) <= eps
            fprintf('[%s] SOC %g → 사용쌍 0 (SIM used: %s)\n', ...
                base_raw, soc_targets(k), strjoin(T.Properties.RowNames(used_idx), ', '));
            continue;
        end

        a_ohm = sum(x.*y) / sum(x.^2);   % Ω
        yhat  = a_ohm*x;
        R2z   = 1 - sum((y - yhat).^2) / max(sum(y.^2), eps);

        R1s_mOhm(k) = a_ohm * 1000;  % mΩ
        R2(k)       = R2z;

        fprintf('[%s] SOC %g: R1=%.3f mΩ, SIM_used=%s, pairs_per_SIM=%s, total n=%d\n', ...
            base_raw, soc_targets(k), R1s_mOhm(k), ...
            strjoin(T.Properties.RowNames(used_idx), ','), ...
            mat2str(pairs_per_SIM{k}.'), n_pairs(k));
    end

    PooledBySOC = table( soc_targets(:), n_pairs(:), R1s_mOhm(:), R2(:), ...
        SIM_used, pairs_per_SIM, ...
        'VariableNames', {'SOC_target','n_pairs','R_1s_mOhm','R2','SIM_used','pairs_per_SIM'} );

    save(fullfile(save_path, sprintf('%s_R1s_pooled_US06_SIM1_SIM9.mat', base_raw)), ...
        'PooledBySOC','soc_targets','soc_tol','sim_pick_names','sim_use_names','-v7.3');

    % 전체 요약 누적
    G_file = [G_file; repmat(string(base_raw), K, 1)];
    G_soc  = [G_soc;  soc_targets(:)];
    G_R1s  = [G_R1s;  R1s_mOhm(:)];
    G_R2   = [G_R2;   R2(:)];
    G_npts = [G_npts; n_pairs(:)];
end

%% 전체 요약 및 피벗 저장 ----------------------------------------------
Summary_all = table(G_file, G_soc, G_R1s, G_R2, G_npts, ...
    'VariableNames', {'file','SOC_target','R_1s_mOhm','R2','n_pairs'});

cellNames = unique(Summary_all.file, 'stable');
socLevels = unique(Summary_all.SOC_target, 'stable');

R1_SOCxCell_mOhm = nan(numel(socLevels), numel(cellNames));
for i = 1:height(Summary_all)
    r = find(socLevels == Summary_all.SOC_target(i));
    c = find(strcmp(cellNames, Summary_all.file(i)));
    R1_SOCxCell_mOhm(r,c) = Summary_all.R_1s_mOhm(i);
end

rowNames1 = arrayfun(@(s) sprintf('SOC %d',s), socLevels, 'uni',false);
varNames1 = matlab.lang.makeValidName(cellNames);
Tbl_SOC_row_Cell_col = array2table(R1_SOCxCell_mOhm, ...
    'RowNames', rowNames1, 'VariableNames', varNames1);

R1_CellxSOC_mOhm = R1_SOCxCell_mOhm.';
rowNames2 = varNames1;
varNames2 = arrayfun(@(s) sprintf('SOC_%d',s), socLevels, 'uni',false);
Tbl_Cell_row_SOC_col = array2table(R1_CellxSOC_mOhm, ...
    'RowNames', rowNames2, 'VariableNames', varNames2);

save(fullfile(save_path, 'R1s_pooled_all_US06_SIM1_SIM9.mat'), ...
    'Summary_all','soc_targets','soc_tol','sim_pick_names', ...
    'R1_SOCxCell_mOhm','R1_CellxSOC_mOhm', ...
    'Tbl_SOC_row_Cell_col','Tbl_Cell_row_SOC_col','-v7.3');

disp('완료: SIM1,SIM9(US06)만 사용한 1초저항 결과 저장 완료.');
