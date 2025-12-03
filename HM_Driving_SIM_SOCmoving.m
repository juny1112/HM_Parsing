% ======================================================================
%  모든 *.mat(‘parsed_data’) → SIM 스텝 요약 테이블
% ----------------------------------------------------------------------
%  • SIM_table : 모든 SIM 스텝 정보
% ======================================================================

clc; clear; close all;

folder_path = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\parsed_data';
save_path   = 'G:\공유 드라이브\BSL_Data4\HNE_agedcell_8_processed\SIM_parsed';

files = dir(fullfile(folder_path,"*.mat"));

for f = 1:numel(files)

    % 1) .mat 로드
    load(fullfile(folder_path,files(f).name),"parsed_data");

    % 2) step_type 클린업
    rawTypes  = string({parsed_data.step_type});
    cleanType = regexprep(rawTypes,"[\s\r\n]+","");

    % 3) SIM 인덱스
    sim_idx = find(strcmpi(cleanType,"SIM"));
    nSIM    = numel(sim_idx);
    if nSIM==0, continue, end

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

        time_raw      = parsed_data(k).time(:);
        timeSIM{s}    = time_raw;

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
            charge_cum   = cumtrapz(t_rel,currentSIM{s});
            charge_total = trapz(t_rel,currentSIM{s});
            if charge_total~=0
                soc_vec = SOC1(s) + (charge_cum/charge_total)*(SOC2(s)-SOC1(s));
            else
                soc_vec = SOC1(s)*ones(size(t_rel));
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

    % 7) 저장(.mat)  → 두 테이블 모두 보관
    [~,base,~] = fileparts(files(f).name);
    save(fullfile(save_path,[base '_SIM.mat']),"SIM_table");
    fprintf("[done] %s → SIM table (%d행) 저장\n", files(f).name, height(SIM_table));
end

fprintf("\n모든 파일 처리 완료!\n");


