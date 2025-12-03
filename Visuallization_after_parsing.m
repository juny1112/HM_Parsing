clc; clear; close all;

% 경로 설정
data_folder = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\parsed_data';
fig_folder  = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\parsed_data\Fig';

if ~exist(fig_folder, 'dir'), mkdir(fig_folder); end

files = dir(fullfile(data_folder, '*.mat'));
N = numel(files);

for k = 1:N
    mat_path = fullfile(files(k).folder, files(k).name);
    S = load(mat_path);  % parsed_data 또는 intrim 중 하나가 있어야 함

    % ---- 데이터 꺼내기 -------------------------------------------------
    if isfield(S, 'parsed_data')
        pd = S.parsed_data;

        % 빈 스텝 제거
        valid = arrayfun(@(s) ~isempty(s.time), pd);
        pd = pd(valid);
        if isempty(pd)
            warning('유효한 스텝이 없습니다. 건너뜀: %s', files(k).name);
            continue;
        end

        % 전 스텝 이어붙이기
        t_all = vertcat(pd.time);      % duration
        v_all = vertcat(pd.voltage);
        i_all = vertcat(pd.current);

        % 스텝 경계(원하면 사용)
        % t_starts = arrayfun(@(s) s.time(1), pd);    % duration
        % xline_pos_sec = seconds(t_starts - t_all(1));

    elseif isfield(S, 'intrim')
        % 예전 형식 호환
        t_all = S.intrim.time;         % duration
        v_all = S.intrim.voltage;
        i_all = S.intrim.current;
        % xline_pos_sec = [];           % (intrim만 있을 때는 생략)
    else
        warning('parsed_data/intrim 모두 없음. 건너뜀: %s', files(k).name);
        continue;
    end

    % 열 벡터화 및 시간축(초)로 변환(시작=0 s)
    t_all = t_all(:); v_all = v_all(:); i_all = i_all(:);
    t_sec = seconds(t_all - t_all(1));

    % ---- 플롯 ----------------------------------------------------------
    fig = figure('Visible','on');  % 배치면 'off' 추천
    yyaxis left
    plot(t_sec, i_all, 'b-', 'LineWidth', 1); hold on
    ylabel('Current [A]')

    yyaxis right
    plot(t_sec, v_all, 'r-', 'LineWidth', 1); hold off
    ylabel('Voltage [V]')

    title(strrep(files(k).name, '_', '\_'), 'Interpreter','tex')
    xlabel('Time [s]')
    box on; grid on;
    set(gcf, 'Position', [300 300 1200 800])

    % (옵션) 스텝 경계선 표시
    % if exist('xline_pos_sec','var') && numel(xline_pos_sec) > 1
    %     hold on; xline(xline_pos_sec(2:end), ':', 'HandleVisibility','off'); hold off
    % end

    % ---- 저장 ----------------------------------------------------------
    [~, base, ~] = fileparts(files(k).name);
    fig_name = fullfile(fig_folder, base + ".fig");
    savefig(fig, fig_name);

    % PNG도 저장하려면 주석 해제
    % print(fig, fullfile(fig_folder, base + ".png"), '-dpng', '-r150');

    %close(fig);
    fprintf("%d/%d 저장 완료: %s\n", k, N, fig_name);
end

fprintf('전체 완료!\n');







