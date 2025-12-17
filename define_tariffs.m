function tariffs = define_tariffs()
    % Define tariff profiles for all locations
    tariffs = struct();

    %% California (PG&E)
    % Summer (June-Sep)
    tariffs.California.summer.weekday.residential = [0.36*ones(1,16), 0.59*ones(1,5), 0.36*ones(1,3)];
    tariffs.California.summer.weekday.commercial = [0.27*ones(1,16), 0.4425*ones(1,5), 0.27*ones(1,3)];
    
    % Winter (Oct-May)
    tariffs.California.winter.weekday.residential = [0.39*ones(1,8), 0.35*ones(1,8), 0.52*ones(1,5), 0.39*ones(1,3)];
    tariffs.California.winter.weekday.commercial = [0.2925*ones(1,8), 0.2625*ones(1,8), 0.39*ones(1,5), 0.2925*ones(1,3)];
    
    % Weekends (assumed same as weekday for California)
    tariffs.California.summer.weekend = tariffs.California.summer.weekday;
    tariffs.California.winter.weekend = tariffs.California.winter.weekday;
    
    %% New York (ConEd)
    % Summer (June-Sep)
    tariffs.NewYork.summer.weekday.residential = [0.025*ones(1,9), 0.35*ones(1,15)]; 
    tariffs.NewYork.summer.weekday.commercial = [0.019*ones(1,9), 0.525*ones(1,15)];
    
    % Winter (Oct-May)
    tariffs.NewYork.winter.weekday.residential = [0.025*ones(1,9), 0.13*ones(1,15)];
    tariffs.NewYork.winter.weekday.commercial = [0.019*ones(1,9), 0.26*ones(1,15)];
    
    % Weekends (all 0.019 for commercial)
    tariffs.NewYork.summer.weekend.residential = tariffs.NewYork.summer.weekday.residential; 
    tariffs.NewYork.summer.weekend.commercial = 0.019*ones(1,24);
    tariffs.NewYork.winter.weekend.residential = tariffs.NewYork.winter.weekday.residential;
    tariffs.NewYork.winter.weekend.commercial = 0.019*ones(1,24);
    
    %% Minnesota
    % Minnesota Power (Residential)
    tariffs.Minnesota.summer.weekday.residential = [0.074*ones(1,5), 0.105*ones(1,10), 0.156*ones(1,5), 0.105*ones(1,3), 0.074];
    tariffs.Minnesota.summer.weekend.residential = [0.074*ones(1,5), 0.105*ones(1,18), 0.074];
    tariffs.Minnesota.winter.weekday.residential = [0.074*ones(1,5), 0.105*ones(1,10), 0.156*ones(1,5), 0.105*ones(1,3), 0.074];
    tariffs.Minnesota.winter.weekend.residential = [0.074*ones(1,5), 0.105*ones(1,18), 0.074];
    
    % XcelEnergy (Commercial) - Year-round
    xc_vals = [0.00810*ones(1,6), 0.02686*ones(1,9), 0.05054*ones(1,5), 0.02686*ones(1,4)];
    tariffs.Minnesota.summer.weekday.commercial = xc_vals;
    tariffs.Minnesota.winter.weekday.commercial = xc_vals;
    tariffs.Minnesota.summer.weekend.commercial = xc_vals;
    tariffs.Minnesota.winter.weekend.commercial = xc_vals;
    % Same for all seasons
    
    %% Texas (PEC)
    % Summer (June-Sep)
    tariffs.Texas.summer.weekday.residential = [0.039905, 0.039905, 0.039905, 0.038387, 0.038387, 0.039905, 0.039905, 0.047026, 0.047026, 0.047026, 0.047026, 0.047026, 0.091961, 0.091961, 0.096305, 0.096305, 0.096305, 0.096305, 0.091961, 0.091961, 0.047026, 0.047026, 0.047026, 0.039905];
    tariffs.Texas.summer.weekday.commercial = [0.071398, 0.071398, 0.071398, 0.069880, 0.069880, 0.071398, 0.071398, 0.078519, 0.078519, 0.078519, 0.078519, 0.078519, 0.123454, 0.123454, 0.127798, 0.127798, 0.127798, 0.127798, 0.123454, 0.123454, 0.078519, 0.078519, 0.078519, 0.071398];
    
    % Winter (Oct-May)
    tariffs.Texas.winter.weekday.residential = [0.046671, 0.046671, 0.044895, 0.044895, 0.046671, 0.061350, 0.061350, 0.061350, 0.052527, 0.052527, 0.052527, 0.052527, 0.052527, 0.052527, 0.052527, 0.052527, 0.061350, 0.061350, 0.061350, 0.052527, 0.052527, 0.052527, 0.052527, 0.046671];
    tariffs.Texas.winter.weekday.commercial = [0.078164, 0.078164, 0.076388, 0.076388, 0.078164, 0.092843, 0.092843, 0.092843, 0.084020, 0.084020, 0.084020, 0.084020, 0.084020, 0.084020, 0.084020, 0.084020, 0.092843, 0.092843, 0.092843, 0.084020, 0.084020, 0.084020, 0.084020, 0.078164];
    
    % Weekends (assumed same as weekday for Texas)
    tariffs.Texas.summer.weekend = tariffs.Texas.summer.weekday;
    tariffs.Texas.winter.weekend = tariffs.Texas.winter.weekday;
    
    % % Feed in tariffs ($/kWh)
    % tariffs.C_PV = struct(...
    %     'California',  struct('residential', 0.075,  'commercial', 0.145), ...
    %     'NewYork',     struct('residential', 0.10,   'commercial', 0.16), ...
    %     'Texas',       struct('residential', 0.09,   'commercial', 0.10), ...
    %     'Minnesota',   struct('residential', 0.13,   'commercial', 0.25) ...
    % );
end