%% =======================================================================
%  HM_1sR_Temp_plot_US06_SOC50
%  - 입력: 각 온도 폴더의 "R1s_pooled_allLoads_600s_summary.mat"
%  - Summary_all에서 US06 & SOC50만 추출
%  - 셀 키 기준으로 온도들 간 공통 셀 매칭
%  - x축: 온도(0/10/20) / xlim [-10 30]
%  - y축: R_1s_mOhm
% =======================================================================
clear; clc; close all;

% (A PATCH) Text/TeX interpreter 에러 원천 차단
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

%% ---------------- 사용자 설정 ----------------
% 각 온도에서 "1초저항_allLoads_600s" 결과 폴더(=save_path)
% ※ 당신 코드에서는 save_path = fullfile(folder_SIM, '1초저항_allLoads_600s')
%    folder_SIM = ...\10degC\이름정렬 이므로, 아래처럼 잡는 게 일반적입니다.

save_paths = [
"G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\0degC\이름정렬\1초저항_allLoads_600s"
"G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\10degC\이름정렬\1초저항_allLoads_600s"
"G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\20degC\이름정렬\1초저항_allLoads_600s"
];

tempsC   = [0 10 20];     % save_paths 순서와 동일해야 함
loadPick = "US06";
socPick  = 50;

% 색상 매핑용 스칼라 (예: QC/40)  ※ 공통 셀 개수와 안 맞으면 자동 보정
capVec = [57.49;57.57;54;52.22;53.45;51.28;57.91;56.51;42.14;57.27;57.18;58.4];
capLabel = "Capacity (QC/40, Ah)";

outDir = "G:\공유 드라이브\BSL_Data4\HNE_RPT_@50,70_251214_9\Driving\SIM_parsed\1sR_temp_plot_US06_SOC50";
filePrefix  = "US06_SOC50_1sR_";
titlePrefix = "R_1s | US06 | SOC50";
doLegend = true;
xlimTemp = [-2 22];

%% ---------------- 로드 & (US06,SOC50) 추출 ----------------
assert(numel(save_paths)==numel(tempsC), "save_paths 개수와 tempsC 개수가 다릅니다.");

D = struct(); % D(ti).T, D(ti).keys, D(ti).R1 (nCell x 1)

for ti = 1:numel(tempsC)
    T = tempsC(ti);

    matPath = fullfile(save_paths(ti), "R1s_pooled_allLoads_600s_summary.mat");
    if ~isfile(matPath)
        error("[%ddegC] summary mat이 없습니다: %s", T, matPath);
    end

    S = load(matPath, "Summary_all");
    if ~isfield(S,"Summary_all") || ~istable(S.Summary_all)
        error("[%ddegC] Summary_all 테이블을 찾지 못했습니다: %s", T, matPath);
    end
    Summary_all = S.Summary_all;

    % 필수 컬럼 체크
    needCols = {'file','load','SOC_target','R_1s_mOhm'};
    if ~all(ismember(needCols, Summary_all.Properties.VariableNames))
        error("[%ddegC] Summary_all에 필요한 컬럼이 없습니다. need=%s", ...
            T, strjoin(needCols, ", "));
    end

    % US06 & SOC50만
    m = strcmp(string(Summary_all.load), string(loadPick)) & (Summary_all.SOC_target == socPick);
    SL = Summary_all(m, :);

    if isempty(SL)
        warning("[%ddegC] US06 & SOC50 데이터가 비었습니다. (해당 온도에서 계산/저장이 되었는지 확인)", T);
        D(ti).T = T;
        D(ti).keys = {};
        D(ti).R1 = [];
        continue;
    end

    % 셀 키 추출 (온도별 파일명/표기 차이 때문에 key로 맞추는 게 안전)
    rawFile = string(SL.file);
    keys = cell(numel(rawFile),1);
    for i = 1:numel(rawFile)
        keys{i} = extractCellKey(rawFile(i));
    end

    R1 = SL.R_1s_mOhm;

    D(ti).T = T;
    D(ti).keys = keys;
    D(ti).R1 = R1(:);
end

%% ---------------- 공통 셀 매칭 (키 기반) ----------------
% 비어있는 온도가 있으면 교집합이 공집합이 될 수 있으므로 먼저 체크
for ti=1:numel(D)
    if isempty(D(ti).keys)
        error("온도 %ddegC에서 (US06,SOC50) 데이터가 비어있어 플롯 불가. 먼저 1초저항 계산을 완료하세요.", D(ti).T);
    end
end

commonKeys = D(1).keys;
for ti = 2:numel(D)
    commonKeys = intersect(commonKeys, D(ti).keys, 'stable');
end

if isempty(commonKeys)
    disp("=== 디버그: 온도별 key 샘플(상위 30개) ===");
    for ti = 1:numel(D)
        fprintf("[%ddegC] keys sample:\n", D(ti).T);
        disp(string(D(ti).keys(1:min(30,numel(D(ti).keys))))');
    end
    error("키 기반 교집합이 비었습니다. extractCellKey 규칙이 file명 패턴과 맞는지 확인 필요.");
end

nCommon = numel(commonKeys);
fprintf(">> 공통 셀 수(키 기준): %d\n", nCommon);

% 공통키 순서로 온도별 R1 재정렬
R1_byT = nan(nCommon, numel(D));  % [cell x temp]
for ti = 1:numel(D)
    [tf, loc] = ismember(commonKeys, D(ti).keys);
    if any(~tf)
        error("내부 오류: 공통키가 특정 온도에서 매칭되지 않았습니다.");
    end
    R1_byT(:, ti) = D(ti).R1(loc);
end

cellNames = commonKeys;

%% ---------------- capVec 길이 보정 ----------------
if numel(capVec) ~= nCommon
    warning("capVec 길이(%d) != 공통 셀 수(%d). min~max 범위로 자동 재생성합니다.", numel(capVec), nCommon);
    if isempty(capVec)
        capVec = linspace(48, 59, nCommon)';
    else
        capVec = linspace(min(capVec), max(capVec), nCommon)';
    end
else
    capVec = capVec(:);
end

%% ---------------- 컬러맵 (기존 스타일 유지) ----------------
Nmap = 256;
anchors = [0.88 0.16 0.24; 0.83 0.70 0.86; 0.16 0.38 0.92];
x  = [0 0.5 1];
xi = linspace(0,1,Nmap)';
cmap = [interp1(x,anchors(:,1),xi,'pchip'), ...
        interp1(x,anchors(:,2),xi,'pchip'), ...
        interp1(x,anchors(:,3),xi,'pchip')];
cmap = min(max(cmap,0),1);
hsvv = rgb2hsv(cmap);
hsvv(:,2) = max(0.35,hsvv(:,2));
hsvv(:,2) = min(1.0,hsvv(:,2)*1.2);
hsvv(:,3) = max(0.75,hsvv(:,3)*0.95);
cmap = hsv2rgb(hsvv);

capMin = min(capVec);
capMax = max(capVec);
mapColor = @(v) cmap( max(1, min(Nmap, 1 + round((v-capMin)/max(capMax-capMin,eps)*(Nmap-1)))), : );

%% ---------------- 플롯 ----------------
if ~exist(outDir,'dir'), mkdir(outDir); end

fig = figure('Color','w','Name','US06 SOC50 1sR vs Temp');
hold on; grid on;

xT = tempsC(:)';

for i = 1:nCommon
    col = mapColor(capVec(i));
    y = R1_byT(i, :);
    plot(xT, y, '-o', ...
        'LineWidth', 1.8, ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'DisplayName', cellNames{i});
end

xlabel('Temperature (°C)');
ylabel('R_1s (mΩ)');
title(titlePrefix, 'Interpreter','none');

% ---- y축은 무조건 0부터 시작 (상한은 데이터 기반 자동) ----
ymax = max(R1_byT(:), [], 'omitnan');
if ~isfinite(ymax) || ymax <= 0
    ymax = 1; % 비정상/전부 NaN/전부 0 이하일 때 안전장치
end
ylim([0, 1.05*ymax]);

xlim(xlimTemp);
xticks(tempsC);

colormap(cmap);
cb = colorbar('Location','eastoutside');
cb.Label.String = capLabel;
clim([capMin capMax]);

if doLegend
    legend('Location','best', 'Interpreter','none');
end

savefig(fig, fullfile(outDir, sprintf('%sR1s_vs_Temp.fig', filePrefix)));
exportgraphics(fig, fullfile(outDir, sprintf('%sR1s_vs_Temp.png', filePrefix)), 'Resolution', 220);

disp("완료: US06, SOC50, 1초저항(R1s) vs 온도 플롯 저장 완료.");

%% ===================== 보조 함수 =====================
function key = extractCellKey(name0)
% Summary_all.file(=base_raw)에서 온도에 무관한 "셀 키"를 뽑는다.
% 1) "01_HNE_..." 처럼 시작 숫자
% 2) "x05_..." 형태
% 3) fallback: 온도 토큰 제거 후 정리

    s = char(string(name0));

    tok = regexp(s, '^(\d+)_', 'tokens', 'once');
    if ~isempty(tok)
        key = tok{1};
        return
    end

    tok = regexp(s, '^x(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        key = tok{1};
        return
    end

    s = regexprep(s, '(0|10|20|30|35|45)degC', '');
    s = regexprep(s, '\s+', '');
    s = regexprep(s, '_{2,}', '_');
    key = s;
end
