clc, clear, close all

% folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_SOC_moving_cutoff_5\raw_data';
% folder_path      = 'G:\공유 드라이브\BSL_Data4\HNE_SOC_moving_cutoff_5\추가';
% save_folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_SOC_moving_cutoff_5_processed\Driving_parsed';

% folder_path      = 'G:\공유 드라이브\BSL_Data4\HNE_Integrated_6';
% save_folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_Integrated_6_processed\parsed data';

% folder_path      = 'G:\공유 드라이브\BSL_Data4\HNE_Integrated_6\order3';
% save_folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_Integrated_6_processed\Test4(order3)\parsed data';

folder_path      = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8\csv파일';
save_folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\parsed_data';

% CSV 파일만 선택
folder = dir(fullfile(folder_path, '*.csv'));
N = length(folder);

for i = 1:N
    % file now assign
    file_name_now = folder(i).name;
    file_path_now = fullfile(folder_path, file_name_now);

    % data now load
    data_now = readtable(file_path_now,"VariableNamingRule","preserve");

    % data now assign
    voltage_now     = data_now.("Voltage(V)");
    current_now     = data_now.("Current(A)");
    time_now        = seconds(data_now.("Total Time"));
    step_time_now   = seconds(data_now.("Time"));
    index_now       = data_now.("DataPoint");
    cycle_index_now = data_now.("Cycle Index");
    % temp_now        = data_now.("T1(℃)");
    step_index_now  = data_now.("Step Index");    % Step Index 기준 파싱
    step_type_now   = data_now.("Step Type");

    % intrim 구조체 저장 (전체)
    intrim.voltage     = voltage_now;
    intrim.current     = current_now;
    intrim.time        = time_now;
    intrim.step_time   = step_time_now;
    intrim.index       = index_now;
    intrim.cycle_index = cycle_index_now;
    % intrim.temp        = temp_now;
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
    disp('Start parsing process')
    parsed_data = struct();
    for j = 1:step_no
        idx_this_step = (intrim.step == j);
        parsed_data(j).voltage      = intrim.voltage(idx_this_step);
        parsed_data(j).current      = intrim.current(idx_this_step);
        parsed_data(j).time         = intrim.time(idx_this_step);
        parsed_data(j).step_time    = intrim.step_time(idx_this_step);
        parsed_data(j).index        = intrim.index(idx_this_step);
        parsed_data(j).cycle_index  = intrim.cycle_index(idx_this_step);
        % parsed_data(j).temp         = intrim.temp(idx_this_step);
        
        parsed_data(j).step_type    = intrim.step_type(idx_this_step);
        parsed_data(j).step_index   = intrim.step_index(idx_this_step);
        parsed_data(j).step_type    = parsed_data(j).step_type(1);
        parsed_data(j).step_index   = parsed_data(j).step_index(1); % 실제 step index
        % parsed_data(j).parsed_step  = j; % 파싱된 순서
    end

    % ---- 저장: 확장자 강제 .mat ----
    [~, base, ~] = fileparts(file_name_now);         % 이름과 확장자 분리
    out_name = base + ".mat";                        % 확장자 .mat로 강제
    out_path = fullfile(save_folder_path, out_name);
    save(out_path, 'parsed_data');                   % 필요시 '-v7.3'

end




