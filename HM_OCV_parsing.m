clc, clear, close all

% OCV 가져오기

parsed_folder = 'C:\Users\junny\OneDrive\문서\GitHub\Parsing\RPT_parsed_data';
save_folder_path = 'C:\Users\junny\OneDrive\문서\GitHub\Parsing\OCV_parsed_data';

if ~exist(save_folder_path, 'dir'), mkdir(save_folder_path); end

parsed_files = dir(fullfile(parsed_folder, '*.mat'));

for k = 1:length(parsed_files)
    file_now = parsed_files(k).name;
    load(fullfile(parsed_folder, file_now), 'parsed_data'); % 'intrim' 불필요

    % 충전 구간 (16번째)
    V_charge = parsed_data(16).voltage;
    N_charge = length(V_charge);
    SOC_charge = linspace(0, 100, N_charge)'; % 0 → 100

    % 방전 구간 (18번째)
    V_discharge = parsed_data(18).voltage;
    N_discharge = length(V_discharge);
    SOC_discharge = linspace(100, 0, N_discharge)'; % 100 → 0

    % 충전(16)
    T_charge = table(SOC_charge, V_charge, 'VariableNames', {'SOC', 'OCV'});
    
    % 방전(18)
    T_discharge = table(SOC_discharge, V_discharge, 'VariableNames', {'SOC', 'OCV'});

    OCV_curve.charge = T_charge;
    OCV_curve.discharge = T_discharge;

    % 저장
    % save(fullfile(save_folder_path, erase(file_now, '.mat') + "_SOC_OCV.mat"), 'OCV_curve');

    % (옵션) 플롯
    figure('Name', file_now)
    plot(OCV_curve.discharge.SOC, OCV_curve.discharge.OCV, 'r.-', 'LineWidth',1.5); hold on
    plot(OCV_curve.charge.SOC, OCV_curve.charge.OCV, 'b.-', 'LineWidth',1.5); hold off
    legend('Discharge', 'Charge')
    xlabel('SOC [%]')
    ylabel('OCV [V]')
    title(['SOC-OCV Curve: ' strrep(file_now, '_', ' ')])
    grid on
end
