clc; clear; close all

%% ====== 경로 설정 ======
folder_path      = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\0degC';
save_folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\parsed_data\0degC';

if ~exist(save_folder_path,'dir'), mkdir(save_folder_path); end

%% ====== 실행 모드 선택 ======
% "all" : 폴더 내 csv 전체 처리
% "one" : 특정 파일 1개만 처리
RUN_MODE = "one"; 

% RUN_MODE="one"일 때 사용할 파일명(확장자 포함)
ONE_FILE_NAME = "HNE_0degC_2C 0.33C_20cyc_18_ch8_0113.csv";   % 예: "SIM1.csv"

%% ====== 처리 대상 파일 리스트 만들기 ======
switch RUN_MODE
    case "all"
        folder = dir(fullfile(folder_path, '*.csv'));
        if isempty(folder)
            error("폴더에 CSV가 없습니다: %s", folder_path);
        end

    case "one"
        file_path_one = fullfile(folder_path, ONE_FILE_NAME);
        if ~isfile(file_path_one)
            error("지정한 파일이 없습니다: %s", file_path_one);
        end
        folder = dir(file_path_one);  % dir 구조체 형태 맞추기 위해
    otherwise
        error('RUN_MODE는 "all" 또는 "one"만 가능합니다.');
end

N = length(folder);

%% ====== 메인 루프 ======
for i = 1:N
    % file now assign
    file_name_now = folder(i).name;
    file_path_now = fullfile(folder_path, file_name_now);

    try
        % data now load
        data_now = readtable(file_path_now,"VariableNamingRule","preserve");

        % data now assign
        voltage_now     = data_now.("Voltage(V)");
        current_now     = data_now.("Current(A)");
        time_now        = seconds(data_now.("Total Time"));
        step_time_now   = seconds(data_now.("Time"));
        index_now       = data_now.("DataPoint");
        cycle_index_now = data_now.("Cycle Index");
        step_index_now  = data_now.("Step Index");    % Step Index 기준 파싱
        step_type_now   = data_now.("Step Type");

        % intrim 구조체 저장 (전체)
        intrim.voltage     = voltage_now;
        intrim.current     = current_now;
        intrim.time        = time_now;
        intrim.step_time   = step_time_now;
        intrim.index       = index_now;
        intrim.cycle_index = cycle_index_now;
        intrim.step_index  = step_index_now;
        intrim.step_type   = step_type_now;

        % Step Index 기준으로 step 부여
        intrim.step = zeros(size(step_index_now));
        intrim.step(1) = 1;
        step_no = 1;
        for n = 2:length(step_index_now)
            if step_index_now(n) == step_index_now(n-1)
                intrim.step(n) = step_no;
            else
                step_no = step_no + 1;
                intrim.step(n) = step_no;
            end
        end

        % Step 별로 parsing
        parsed_data = struct();
        for j = 1:step_no
            idx_this_step = (intrim.step == j);
            parsed_data(j).voltage      = intrim.voltage(idx_this_step);
            parsed_data(j).current      = intrim.current(idx_this_step);
            parsed_data(j).time         = intrim.time(idx_this_step);
            parsed_data(j).step_time    = intrim.step_time(idx_this_step);
            parsed_data(j).index        = intrim.index(idx_this_step);
            parsed_data(j).cycle_index  = intrim.cycle_index(idx_this_step);

            parsed_data(j).step_type    = intrim.step_type(idx_this_step);
            parsed_data(j).step_index   = intrim.step_index(idx_this_step);

            parsed_data(j).step_type  = parsed_data(j).step_type(1);
            parsed_data(j).step_index = parsed_data(j).step_index(1); % 실제 step index
        end

        % ---- 저장: 확장자 강제 .mat ----
        [~, base, ~] = fileparts(file_name_now);
        out_name = base + ".mat";
        out_path = fullfile(save_folder_path, out_name);
        save(out_path, 'parsed_data');   % 필요시 '-v7.3'

        fprintf("[%d/%d] Saved: %s\n", i, N, out_name);

    catch ME
        warning("[%d/%d] FAIL: %s | %s", i, N, file_name_now, ME.message);
        continue;
    end
end

disp("Done.")
