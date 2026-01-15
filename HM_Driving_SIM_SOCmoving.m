% ======================================================================
%  (선택 모드) 모든 *.mat(‘parsed_data’) 또는 단일 파일 → SIM 스텝 요약 테이블
% ----------------------------------------------------------------------
%  • 입력: parsed_data가 들어있는 *.mat
%  • 출력: *_SIM.mat (SIM_table)
% ======================================================================

clc; clear; close all;

%% ====== 경로 설정 ======
folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\parsed_data\0degC';
save_path   = 'G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\0degC';

if ~exist(save_path,'dir'), mkdir(save_path); end

%% ====== 실행 모드 선택 ======
% "all" : 폴더 내 *.mat 전체 처리
% "one" : 특정 *.mat 1개만 처리
RUN_MODE = "one";        
ONE_FILE_NAME = "HNE_0degC_2C 0.33C_20cyc_18_ch8_0113.mat"; 

%% ====== 처리 대상 파일 리스트 만들기 ======
switch RUN_MODE
    case "all"
        files = dir(fullfile(folder_path,"*.mat"));
        if isempty(files)
            error("폴더에 .mat 파일이 없습니다: %s", folder_path);
        end

    case "one"
        file_one = fullfile(folder_path, ONE_FILE_NAME);
        if ~isfile(file_one)
            error("지정한 파일이 없습니다: %s", file_one);
        end
        files = dir(file_one);  % dir 구조 유지

    otherwise
        error('RUN_MODE는 "all" 또는 "one"만 가능합니다.');
end

fprintf("Target files: %d (%s)\n\n", numel(files), RUN_MODE);

%% ====== 메인 루프 ======
for f = 1:numel(files)

    % 1) .mat 로드 (parsed_data 존재 확인)
    mat_path = fullfile(folder_path,files(f).name);

    S = load(mat_path); % 일단 전체 로드 후 필드 체크
    if ~isfield(S,"parsed_data")
        warning("[skip] %s : parsed_data 변수가 없음", files(f).name);
        continue
    end
    parsed_data = S.parsed_data;

    % 2) step_type 클린업
    rawTypes  = string({parsed_data.step_type});
    cleanType = regexprep(rawTypes,"[\s\r\n]+","");

    % 3) SIM 인덱스
    sim_idx = find(strcmpi(cleanType,"SIM"));
    nSIM    = numel(sim_idx);
    if nSIM==0
        fprintf("[skip] %s : SIM step 없음\n", files(f).name);
        continue
    end

    % 4) 결과 변수
    voltageSIM = cell(nSIM,1);
    currentSIM = cell(nSIM,1);
    timeSIM    = cell(nSIM,1);
    SOC_vecSIM = cell(nSIM,1);
    OCV_vecSIM = cell(nSIM,1);
    OCV1 = nan(nSIM,1);  OCV2 = nan(nSIM,1);
    SOC1 = nan(nSIM,1);  SOC2 = nan(nSIM,1);

    rest_mask = strcmpi(cleanType,"Rest");
    idxVec    = 1:numel(parsed_data);

    % 5) SIM별 계산
    for s = 1:nSIM
        k = sim_idx(s);

        voltageSIM{s} = parsed_data(k).voltage(:);
        currentSIM{s} = parsed_data(k).current(:);

        time_raw   = parsed_data(k).time(:);
        timeSIM{s} = time_raw;

        % Rest 전·후
        prevRest = find(rest_mask & (idxVec < k),1,"last");
        if ~isempty(prevRest)
            OCV1(s) = parsed_data(prevRest).voltage(end);
            SOC1(s) = SOC_interp(OCV1(s));
        end
        nextRest = find(rest_mask & (idxVec > k),1,"first");
        if ~isempty(nextRest)
            OCV2(s) = parsed_data(nextRest).voltage(end);
            SOC2(s) = SOC_interp(OCV2(s));
        end

        % 쿨롱카운팅
        if ~isnan(SOC1(s)) && ~isnan(SOC2(s))
            t_rel        = time_raw - time_raw(1);
            charge_cum   = cumtrapz(t_rel, currentSIM{s});
            charge_total = trapz(t_rel, currentSIM{s});

            if charge_total ~= 0
                soc_vec = SOC1(s) + (charge_cum/charge_total) * (SOC2(s)-SOC1(s));
            else
                soc_vec = SOC1(s) * ones(size(t_rel));
            end
        else
            soc_vec = nan(size(time_raw));
        end

        SOC_vecSIM{s} = soc_vec;
        OCV_vecSIM{s} = OCV_interp(soc_vec);
    end

    % 6) 전체 SIM 테이블
    rowNames = cellstr("SIM"+string((1:nSIM).'));
    varNames = {'voltage','current','time', ...
                'OCV1','OCV2','SOC1','SOC2','SOC_vec','OCV_vec'};

    SIM_table = table( ...
        voltageSIM, currentSIM, timeSIM, ...
        OCV1, OCV2, SOC1, SOC2, ...
        SOC_vecSIM, OCV_vecSIM, ...
        'VariableNames',varNames, ...
        'RowNames',rowNames );

    % 7) 저장(.mat)
    [~,base,~] = fileparts(files(f).name);
    out_file = fullfile(save_path, [base '_SIM.mat']);
    save(out_file, "SIM_table");

    fprintf("[done] %s → SIM table (%d행) 저장\n", files(f).name, height(SIM_table));
end

fprintf("\n모든 파일 처리 완료!\n");
