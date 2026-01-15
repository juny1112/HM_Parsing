% ======================================================================
%  SIM_table: time-vs-current + '주행부하_cutoff' 템플릿 오버레이 (4x8)
%  - 파란 실선: 실제 전류 (SIM)
%  - 빨간 점선: 인가 전류 (프로토콜)
%  - 범례: 우상단, 그림당 1개
%  - 저장: ...\SIM_parsed\전류인가
% ======================================================================

clc; clear; close all;

% ── 경로 설정 -----------------------------------------------------------
sim_path      = 'G:\공유 드라이브\BSL_Data4\HNE_fresh_integrated_7_Drivingprocessed\SIM_parsed';
save_dir      = fullfile(sim_path, '전류인가');   % 요청 폴더
protocol_path = 'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_cutoff';

if ~exist(save_dir,'dir'); mkdir(save_dir); end

sim_files = dir(fullfile(sim_path, '*_SIM.mat'));
if isempty(sim_files)
    error("SIM 파일(*_SIM.mat)을 찾을 수 없습니다: %s", sim_path);
end

% ── 템플릿(엑셀) 로드 ---------------------------------------------------
tmpl_keys = {'US06','UDDS','HWFET','WLTP','CITY1','CITY2','HW1','HW2'};
xlsx_list = dir(fullfile(protocol_path, '*.xls*'));
if isempty(xlsx_list)
    error("엑셀 템플릿 파일을 찾지 못했습니다: %s", protocol_path);
end

find_cols = @(T) deal( ...
    find(ismember(lower(string(T.Properties.VariableNames)), ...
         lower(["time","t","sec","seconds"])) ,1,"first"), ...
    find(ismember(lower(string(T.Properties.VariableNames)), ...
         lower(["current","i","amp","a"])) ,1,"first") );
contains_ci = @(s,kw) contains(lower(s), lower(kw));

templates = struct('name',{},'t',{},'i',{});
for k = 1 %1:numel(tmpl_keys)
    key = tmpl_keys{k};
    match = [];
    for j = 1:numel(xlsx_list)
        if contains_ci(xlsx_list(j).name, key)
            match = xlsx_list(j); break
        end
    end
    if isempty(match)
        warning("키워드 '%s'에 해당하는 엑셀을 찾지 못했습니다. 건너뜁니다.", key);
        continue
    end

    T = readtable(fullfile(protocol_path, match.name));
    [t_idx, i_idx] = find_cols(T);
    if isempty(t_idx) || isempty(i_idx), t_idx = 1; i_idx = 2; end

    tt = T{:, t_idx};
    ii = T{:, i_idx};

    if isdatetime(tt),    tt = seconds(tt - tt(1));
    elseif isduration(tt),tt = seconds(tt - tt(1));
    else,                 tt = tt - tt(1);
    end

    templates(end+1).name = key; %#ok<SAGROW>
    templates(end).t = tt(:);
    templates(end).i = ii(:);
end
if numel(templates) < 8
    warning('템플릿이 8종 미만입니다(%d). 있는 것만 사용합니다.', numel(templates));
end

% ── 서브플롯/스케일 옵션 -----------------------------------------------
nrows = 4; ncols = 8; perFig = nrows * ncols;
time_unit = 's';          % 's' | 'min' | 'h'
amp_match = true;         % 템플릿 전류 진폭을 SIM과 대략 맞출지

to_rel_seconds = @(t) ...
    ( isdatetime(t) * seconds(t - t(1)) ) + ...
    ( isduration(t) * seconds(t - t(1)) ) + ...
    ( isnumeric(t)  * (t - t(1)) );

% ── 메인 루프 -----------------------------------------------------------
for f = 1:numel(sim_files)
    S = load(fullfile(sim_path, sim_files(f).name), 'SIM_table');
    if ~isfield(S, 'SIM_table')
        warning("'%s'에 SIM_table 변수가 없습니다. 건너뜁니다.", sim_files(f).name);
        continue;
    end
    T = S.SIM_table;
    nSIM = height(T);
    if nSIM == 0
        warning("'%s'의 SIM 개수가 0입니다. 건너뜁니다.", sim_files(f).name);
        continue;
    end

    nPage = ceil(nSIM / perFig);
    for p = 1:nPage
        idx_start = (p-1)*perFig + 1;
        idx_end   = min(p*perFig, nSIM);

        hfig = figure('Color','w','Position',[60 60 1800 950]);
        legend_added = false;   % 그림당 1회만 범례

        for k = idx_start:idx_end
            axIdx = k - idx_start + 1;
            ax = subplot(nrows, ncols, axIdx);

            % --- SIM 데이터 ---
            ts = T.time{k};
            Is = T.current{k};
            if isempty(ts) || isempty(Is)
                title(sprintf('SIM %d (empty)', k), 'FontSize', 9, 'Interpreter','none');
                axis off; continue;
            end

            tsec = to_rel_seconds(ts);
            switch lower(time_unit)
                case 's',   tx = tsec;      xlab = 'Time (s)';
                case 'min', tx = tsec/60;   xlab = 'Time (min)';
                case 'h',   tx = tsec/3600; xlab = 'Time (h)';
                otherwise,  tx = tsec;      xlab = 'Time (s)';
            end

            hReal = plot(ax, tx, Is, 'b-', 'LineWidth', 0.9); hold(ax, 'on');
            grid(ax, 'on'); box(ax, 'on');
            xlim(ax, [tx(1), tx(end)]);

            % --- 템플릿 선택(8종 반복) ---
            tmpl_idx = mod(k-1, numel(templates)) + 1;
            name_tag = 'N/A'; hProt = [];

            if tmpl_idx <= numel(templates)
                tt  = templates(tmpl_idx).t;
                ii  = templates(tmpl_idx).i;
                name_tag = templates(tmpl_idx).name;

                if numel(tt) >= 2 && any(isfinite(ii))
                    Ts = tsec(end) - tsec(1);
                    tte = tt(end);
                    if tte > 0
                        scale_t = Ts / tte;
                        ii_on_sim = interp1(tt * scale_t, ii, tsec, 'linear', 'extrap');

                        if amp_match
                            num = prctile(abs(Is), 95);
                            den = max(prctile(abs(ii_on_sim), 95), eps);
                            scale_a = num / den;
                        else
                            scale_a = 1.0;
                        end

                        switch lower(time_unit)
                            case 's',   tx_tmpl = tsec;
                            case 'min', tx_tmpl = tsec/60;
                            case 'h',   tx_tmpl = tsec/3600;
                        end

                        hProt = plot(ax, tx_tmpl, scale_a * ii_on_sim, 'r--', 'LineWidth', 1.0);
                    end
                end
            end

            title(ax, sprintf('SIM %d — %s', k, name_tag), 'FontSize', 9, 'Interpreter','none');
            xlabel(ax, xlab); ylabel(ax, 'Current (A)');

            % --- 범례: 첫 번째 서브플롯에만 생성, 우상단 배치 -------------
            if ~legend_added
                if ~isempty(hProt) && isvalid(hProt)
                    legend(ax, [hReal, hProt], {'실제 전류','인가 전류'}, ...
                           'Location','northeast', 'Interpreter','none');
                else
                    legend(ax, hReal, {'실제 전류'}, ...
                           'Location','northeast', 'Interpreter','none');
                end
                legend_added = true;
            end
        end

        sgtitle(sprintf('%s — SIM current vs. time with templates (page %d/%d)', ...
            sim_files(f).name, p, nPage), 'Interpreter','none','FontWeight','bold');

        % 저장
        [~, base, ~] = fileparts(sim_files(f).name);
        out_png = fullfile(save_dir, sprintf('%s_SIM_I_overlay_p%02d.png', base, p));
        try
            exportgraphics(hfig, out_png, 'Resolution', 220);
        catch
            saveas(hfig, out_png);
        end
    end
end
