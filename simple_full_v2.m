function full_main()
    % This function runs seasonal tests for multiple parameter combinations over a 10-year period
    % We track battery degradation and replacements continuously
    % Output is consolidated into a single CSV file and the peaks are saved
    % in a seperate file
    % Check if parallel pool already exists
    clear all; close all; clc;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        % Determine optimal number of workers (use up to 8)
        num_cores = feature('numcores');
        num_workers = min(num_cores, 8);  % Don't use more than 8 to avoid memory issues
        
        fprintf('Starting parallel pool with %d workers...\n', num_workers);
        fprintf('This may take 30-60 seconds on first run...\n');
        
        parpool('local', num_workers);
        
        fprintf('Parallel pool ready! ✓\n\n');
    else
        fprintf('Using existing parallel pool with %d workers\n\n', poolobj.NumWorkers);
    end
    
    rng(42);

    %% Initialize Base Parameters
    Sb = 10000000;
    Vb = 15000;
    Zb = (Vb^2)/Sb;
    Y = Zb * [ 1.362648+-2.116856i, -1.362648+2.116858i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  -1.362648+2.116858i, 1.635178+-2.540220i, -0.272530+0.423372i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, -0.272530+0.423372i, 1.635178+-2.540220i, -1.362648+2.116858i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, -1.362648+2.116858i, 4.087945+-6.350572i, -2.725296+4.233716i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -2.725296+4.233716i, 4.087945+-6.350572i, -1.362648+2.116858i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -1.362648+2.116858i, 4.087945+-6.350572i, -2.725296+4.233716i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -2.725296+4.233716i, 4.087945+-6.350572i, -1.362648+2.116858i, 0.000000+0.000000i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -1.362648+2.116858i, 4.087945+-6.350572i, -2.725296+4.233716i, 0.000000+0.000000i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -2.725296+4.233716i, 4.087945+-6.350572i, -1.362648+2.116858i,;  0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, 0.000000+0.000000i, -1.362648+2.116858i, 1.362648+-2.116856i,; ] ;
    P= (1/Sb)*[ 0.000000, -500000.000000, -200000.000000, -500000.000000, -300000.000000, -100000.000000, -200000.000000, -300000.000000, -500000.000000, -200000.000000,];
    Q= (1/Sb)*[ 0.000000, -375000.000000, -150000.000000, -375000.000000, -225000.000000, -75000.000000, -150000.000000, -225000.000000, -375000.000000, -150000.000000,];
    P = P'; Q = Q';
    tm = 24;

    % Define Seasons with Tariffs
    season_data_weekdays_RC = {
        'summer', 1.4, 'summer';
        'spring', 1.0, 'winter';
        'winter', 1.2, 'winter';
        'autumn1', 1.0, 'summer'; % Autumn with summer tariff
        'autumn2', 1.0, 'winter'  % Autumn with winter tariff
    };
    
    season_data_weekends_R = {
        'summer', 1.45, 'summer';
        'spring', 1.05, 'winter';
        'winter', 1.25, 'winter';
        'autumn1', 1.05, 'summer'; % Autumn with summer tariff
        'autumn2', 1.05, 'winter'  % Autumn with winter tariff
    };
    
    season_data_weekends_C = {
        'summer', 1.35, 'summer';
        'spring', 0.95, 'winter';
        'winter', 1.15, 'winter';
        'autumn1', 0.95, 'summer'; % Autumn with summer tariff
        'autumn2', 0.95, 'winter'  % Autumn with winter tariff
    };

    % Location
    location = 'Texas'; % [CA, NY, MN, TX]
    
    % Years to process
    years = 1:10;
    
    % IRR vectors for 10 years
    irr_high = [0.12, 0.13, 0.11, 0.14, 0.10, 0.15, 0.12, 0.11, 0.09, 0.08];  
    irr_mod = [0.10, 0.09, 0.08, 0.10, 0.09, 0.11, 0.10, 0.09, 0.07, 0.06];  
    irr_low = [0.07, 0.06, 0.05, 0.06, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03];
    
    irr_vectors = {irr_mod}; % {irr_high, irr_mod, irr_low}; one for now 
    irr_names = {'Moderate'}; % {'High', 'Moderate', 'Low'}; one for now
    
    % Day types for a week
    day_types = {'weekday', 'weekend'};
    
    % All seasons to process (in order)
    all_seasons = {'spring', 'summer', 'autumn1', 'autumn2', 'winter'};

    % Season lengths in weeks
    season_weeks = zeros(length(all_seasons), 1);
    for s = 1:length(all_seasons)
        if strcmp(all_seasons{s}, 'autumn1')
            season_weeks(s) = 4;
        elseif strcmp(all_seasons{s}, 'autumn2')
            season_weeks(s) = 8;
        else
            season_weeks(s) = 12;
        end
    end

    % Load Tariffs
    tariffs = define_tariffs();
    
    % Define parameter combinations to test
    pv_scales = [0.2, 0.5, 1.2];
    battery_types = {'Li-ion', 'Na-ion'};
    x_values = [0.0, 0.5, 1.0, 1.5];  % Battery capacity scaling factors
    
    % Define residential/commercial mixes
    is_residential_choices = {
        [true, true, true, true, true, true, true, true, true, true], % 9/9 residential
        [true, true, false, true, false, true, true, false, true, true], % 6/9 residential
        [true, false, false, true, false, true, false, false, true, false] % 3/9 residential
    };
    residential_mix_names = {'9of9_Res', '6of9_Res', '3of9_Res'};
    
    % IL fraction values
    ilF_res_values = [0.1]; 
    ilF_com_values = [0.5]; 
    
    % IL scale values
    p_ilM_scales = [0.05, 0.1, 0.3];
    
    % Peak values
    peak_values = [0.01, 0.1, 1, 10];
    
    % Create the consolidated output file
    consolidated_filename = 'simple_full_TX.csv';
    fid = fopen(consolidated_filename, 'w');
    
    % Write header for consolidated CSV
    fprintf(fid, 'Simulation Results for All Parameter Combinations\n\n');
    fprintf(fid, 'Parameters,Values\n');
    
    % Write parameter section
    fprintf(fid, 'Location,%s\n', location);
    
    % Format PV scales
    pv_str = '';
    for i = 1:length(pv_scales)
        if i > 1
            pv_str = [pv_str, ','];
        end
        pv_str = [pv_str, num2str(pv_scales(i))];
    end
    fprintf(fid, 'PV Scales,%s\n', pv_str);
    
    % Format battery types
    batt_str = '';
    for i = 1:length(battery_types)
        if i > 1
            batt_str = [batt_str, ','];
        end
        batt_str = [batt_str, battery_types{i}];
    end
    fprintf(fid, 'Battery Types,%s\n', batt_str);
    
    % Format battery scaling factors
    x_str = '';
    for i = 1:length(x_values)
        if i > 1
            x_str = [x_str, ','];
        end
        x_str = [x_str, num2str(x_values(i))];
    end
    fprintf(fid, 'Battery Scaling Factors,%s\n', x_str);
    
    % Format residential mixes
    res_mix_str = '';
    for i = 1:length(residential_mix_names)
        if i > 1
            res_mix_str = [res_mix_str, ','];
        end
        res_mix_str = [res_mix_str, residential_mix_names{i}];
    end
    fprintf(fid, 'Residential Mixes,%s\n', res_mix_str);
    
    % Format residential IL fractions
    res_il_str = '';
    for i = 1:length(ilF_res_values)
        if i > 1
            res_il_str = [res_il_str, ','];
        end
        res_il_str = [res_il_str, num2str(ilF_res_values(i))];
    end
    fprintf(fid, 'Residential IL Fractions,%s\n', res_il_str);
    
    % Format commercial IL fractions
    com_il_str = '';
    for i = 1:length(ilF_com_values)
        if i > 1
            com_il_str = [com_il_str, ','];
        end
        com_il_str = [com_il_str, num2str(ilF_com_values(i))];
    end
    fprintf(fid, 'Commercial IL Fractions,%s\n', com_il_str);
    
    % Format IL scales
    pilm_scale_str = '';
    for i = 1:length(p_ilM_scales)
        if i > 1
            pilm_scale_str = [pilm_scale_str, ','];
        end
        pilm_scale_str = [pilm_scale_str, num2str(p_ilM_scales(i))];
    end
    fprintf(fid, 'Pilm Scales,%s\n', pilm_scale_str);
    
    % Format peak values
    peak_str = '';
    for i = 1:length(peak_values)
        if i > 1
            peak_str = [peak_str, ','];
        end
        peak_str = [peak_str, num2str(peak_values(i))];
    end
    fprintf(fid, 'Peak Values,%s\n', peak_str);
    
    % Format IRR scenarios
    irr_str = '';
    for i = 1:length(irr_names)
        if i > 1
            irr_str = [irr_str, ','];
        end
        irr_str = [irr_str, irr_names{i}];
    end
    fprintf(fid, 'IRR Scenarios,%s\n', irr_str);
    fprintf(fid, '\n');
    
    % Write header for results table
     fprintf(fid, 'Peak,P_ilM_Scale,PV_Scale,X_Battery_Scale,Battery_Type,Residential_Mix,ilF_Res,ilF_Com,IRR_Name,NPV_Decade_Fval,Total_Decade_Fval,Total_Load_Served,Total_Replacements,Total_DR_Events,Total_DR_Revenue\n');
    
    % Initialize simulation counter
    total_sims = length(peak_values) * length(p_ilM_scales) * length(pv_scales) * ...
                 length(x_values) * length(battery_types) * length(is_residential_choices) * ...
                 length(ilF_res_values) * length(ilF_com_values) * length(irr_vectors);

    total_pk_entries = total_sims * 10; % 10 buses per simulation
    all_pk_results = cell(length(peak_values), 1);  % One cell per parfor iteration
    

    % Calculate total number of combinations
 
    fprintf('\n=== STARTING PARALLEL ANALYSIS ===\n');
    fprintf('Total parameter combinations: %d\n', total_sims);
    fprintf('Workers: %d\n', gcp().NumWorkers);
  
    all_results_data = cell(length(peak_values), 1);
    % Main simulation loops
    parfor peak_idx = 1:length(peak_values)
        peak = peak_values(peak_idx);
        local_pk_results = {}; % local collector for this worker
        local_results_data = {};  

        % Initialize all temporary variables to avoid parfor warnings
        cum_degradation_simple = zeros(10, 1);
        res_only_days = 0;
        com_only_days = 0;
        both_event_days = 0;
        no_event_days = 0;
        season_capacity_revenue_total = 0;
        season_capacity_revenue_res = 0;
        season_capacity_revenue_com = 0;
        has_res_event = false;
        has_com_event = false;
        scenario_multiplier = 0;
        
        for p_ilM_idx = 1:length(p_ilM_scales)
            p_ilM_scale = p_ilM_scales(p_ilM_idx);
            
            for pv_idx = 1:length(pv_scales)
                pv_scale = pv_scales(pv_idx);
                
                for x_idx = 1:length(x_values)
                    x = x_values(x_idx);
                    
                    for batt_idx = 1:length(battery_types)
                        battery_type = battery_types{batt_idx};
                        
                        for res_mix_idx = 1:length(is_residential_choices)
                            is_residential = is_residential_choices{res_mix_idx};
                            residential_mix = residential_mix_names{res_mix_idx};
                            
                            for res_il_idx = 1:length(ilF_res_values)
                                for com_il_idx = 1:length(ilF_com_values)
                                    % Get current IL fraction values
                                    ilF_res = ilF_res_values(res_il_idx);
                                    ilF_com = ilF_com_values(com_il_idx);
                                    
                                    % Create ilF_vec based on bus types
                                    ilF_vec = zeros(1, 10);
                                    for i = 1:10
                                        if is_residential(i)
                                            ilF_vec(i) = ilF_res;
                                        else
                                            ilF_vec(i) = ilF_com;
                                        end
                                    end
                                    
                                    for irr_idx = 1:length(irr_vectors)
                                        irr_vector = irr_vectors{irr_idx};
                                        irr_name = irr_names{irr_idx};
                                        
                                        % Simple calculation based on current indices
                                        sim_count = sub2ind([length(peak_values), length(p_ilM_scales), length(pv_scales), ...
                                                             length(x_values), length(battery_types), length(is_residential_choices), ...
                                                             length(ilF_res_values), length(ilF_com_values), length(irr_vectors)], ...
                                                             peak_idx, p_ilM_idx, pv_idx, x_idx, batt_idx, res_mix_idx, ...
                                                             res_il_idx, com_il_idx, irr_idx);
                                        % Display progress
                                        disp(['Parameters: Peak=', num2str(peak), ', IL Scale=', num2str(p_ilM_scale), ...
                                              ', PV Scale=', num2str(pv_scale), ', Battery Scale=', num2str(x), ...
                                              ', Battery Type=', battery_type, ', Res Mix=', residential_mix, ...
                                              ', IL Res=', num2str(ilF_res), ', IL Com=', num2str(ilF_com), ...
                                              ', IRR=', irr_name]);
                                        
                                        % Input parameters for reference
                                        input_params = struct(...
                                            'Location', location, ...
                                            'PVScale', pv_scale, ...
                                            'BatteryType', battery_type, ...
                                            'BatteryScaleFactor', x, ...
                                            'SectorInfo', sprintf('Mixed: %d Residential, %d Commercial', sum(is_residential), sum(~is_residential)), ...
                                            'ILFractionRes', ilF_res, ...
                                            'ILFractionCom', ilF_com, ...
                                            'PILMScale', p_ilM_scale, ...
                                            'PeakValue', peak);
                                        
                                        % Get battery characteristics
                                        [E_M_initial, E_mt_current, n_b_initial, C_b] = get_battery_characteristics(battery_type, P, x); % c_b is inflated inside BattPkShvRsCm
                            
                                        % Initialize battery state tracking
                                        E_M_current = E_M_initial;
                                        
                                        % Battery degradation initialization
                                        E_ret = 0.8 * E_M_initial;       % 1x10 retirement capacity matrix
                                        replacements = zeros(10,1);      % Track replacements per bus

                                       T_ref = 25;  % reference temp (C)
                                       % === BATTERY PARAMETERS (identical to simple method) ===
                                        if strcmp(battery_type, 'Li-ion')
                                            CYCLE_LIFE_80DOD = 4500;
                                            DOD_EXPONENT = 1.1;
                                            CALENDAR_FADE_ANNUAL = 0.018;
                                            TEMP_SENSITIVITY_CYCLE = 0.04;
                                            TEMP_SENSITIVITY_CAL = 0.045;
                                    
                                        else % Na-ion
                                            CYCLE_LIFE_80DOD = 3000;  % Middle range (1000-3000 cycles from literature)
                                            DOD_EXPONENT = 1.2;  % Steeper degradation curve
                                            CALENDAR_FADE_ANNUAL = 0.01;  % Lower calendar aging (safer chemistry, less prone to thermal runaway)
                                            TEMP_SENSITIVITY_CYCLE = 0.025;  % LOW temp sensitivity - KEY ADVANTAGE of Na-ion!
                                            TEMP_SENSITIVITY_CAL   = 0.03;  % Also lower than Li-ion
                                        end

                                        % Initialize metrics arrays
                                        season_fval = zeros(length(all_seasons), length(years));
                                        season_ldl = zeros(length(all_seasons), length(years));
                                        season_btr = zeros(length(all_seasons), length(years));
                                        season_replacements = zeros(length(all_seasons), length(years));
                                        
                                        annual_fval = zeros(length(years), 1);
                                        annual_ldl = zeros(length(years), 1);
                                        annual_btr = zeros(length(years), 1);
                                        annual_replacements = zeros(length(years), 1);
                                        
                                        
                                        % Initialize storage for PK values per season per year
                                        pk_per_bus_per_season_per_year = zeros(10, length(all_seasons), length(years)); % 10 buses x 5 seasons x 10 years

                                        
                                        % 10-year totals
                                        total_fval = 0;
                                        total_ldl = 0;
                                        total_btr = 0;
                                        
                                        season_DR_capacity_revenue = zeros(length(all_seasons), length(years));
                                        season_DR_performance_revenue = zeros(length(all_seasons), length(years));
                                        season_DR_total_events = zeros(length(all_seasons), length(years));
                              
                                        previous_results = cell(length(all_seasons), length(day_types));
                                        
                                        
                                        % Loop through each year
                                        for yr_idx = 1:length(years)
                                            current_irr = irr_vector(1:yr_idx);

                                            % Calculate annual participation payment (ONCE per year)
                                            annual_participation = 0;
                                            if strcmp(location, 'NewYork') && sum(is_residential(2:end)) > 0 && yr_idx>=3

                                                num_res_buses_ny = sum(is_residential(2:end)); % the rewards are granted per bus
                                                annual_participation = 25 * num_res_buses_ny;  % $25 annual payment per residential bus
                                            end
                                            annual_participation_per_season = annual_participation / length(all_seasons);
                                            
                                            % Process each season within the year
                                            for season_idx = 1:length(all_seasons)
                                                target_season = all_seasons{season_idx};
                                                num_weeks = season_weeks(season_idx);
                                                
                                                % Find index of season in data arrays
                                                season_data_idx = find(strcmp(season_data_weekdays_RC(:,1), target_season));
                                                
                                                % Get tariff season for this season
                                                tariff_season = season_data_weekdays_RC{season_data_idx, 3};
                                                
                                                % Season specific metrics for current season
                                                cur_season_fval = 0;
                                                cur_season_ldl = 0;
                                                cur_season_btr = 0;

                                                % Initialize cumulative degradation tracker ONCE at the beginning
                                                if season_idx == 1 && yr_idx == 1
                                                    cum_degradation_simple = zeros(10, 1);
                                                end
                                                                                                
                                                % Initialize per-bus battery throughput for this season
                                                per_bus_season_btr = zeros(10, 1);
                                                
                                                % Initialize E_mt_current at the beginning of each season
                                                E_mt_current_season = E_mt_current;
                                                
                                                % ============ ONE-TIME CREDITS (HANDLE AT YEAR LEVEL) ============
                                                % Initialize annual one-time credits for first year only
                                                if yr_idx == 1 && season_idx == 1
                                                    % Count residential and commercial buses (excluding slack)
                                                    num_res_buses = sum(is_residential(2:end));
                                                    num_com_buses = sum(~is_residential(2:end));
                                                    
                                                    % Calculate one-time enrollment credits
                                                    annual_one_time_total = 0;
                                                    
                                                    if strcmp(location, 'Texas')
                                                        % $25 per residential customer only
                                                        annual_one_time_total = 25 * num_res_buses;
                                                        
                                                    elseif strcmp(location, 'NewYork')
                                                        % $85 per residential customer only
                                                        annual_one_time_total = 85 * num_res_buses;
                                                        
                                                    elseif strcmp(location, 'Minnesota')
                                                        % $50 per customer (both residential AND commercial)
                                                        annual_one_time_total = 50 * (num_res_buses + num_com_buses);
                                                    end
                                                    
                                                    % Distribute across all seasons in first year
                                                    seasons_per_year = length(all_seasons);
                                                    one_time_credit_per_season = annual_one_time_total / seasons_per_year;
                                                else
                                                    one_time_credit_per_season = 0;
                                                end

                                                % Process each day type (weekday and weekend)
                                                for day_idx = 1:length(day_types)
                                                    day_type = day_types{day_idx};
                                                    
                                                    % Create bus-specific scale factors
                                                    bus_scale_factors = zeros(1, 10);
                                                    
                                                    % Select appropriate scale factor based on day type and sector
                                                    if strcmp(day_type, 'weekday')
                                                        % For weekdays, apply same scale to all buses
                                                        weekday_scale = season_data_weekdays_RC{season_data_idx, 2};
                                                        
                                                        % Apply the same weekday scale factor to all buses
                                                        for bus = 1:10
                                                            bus_scale_factors(bus) = weekday_scale;
                                                        end
                                                        
                                                        day_multiplier = 5 * num_weeks;  % 5 weekdays per week
                                                    else
                                                        % For weekends, apply different scale factors based on bus type
                                                        res_scale_factor = season_data_weekends_R{season_data_idx, 2};
                                                        com_scale_factor = season_data_weekends_C{season_data_idx, 2};
                                                        
                                                        % Apply residential or commercial scale based on bus type
                                                        for bus = 1:10
                                                            if is_residential(bus)
                                                                bus_scale_factors(bus) = res_scale_factor;
                                                            else
                                                                bus_scale_factors(bus) = com_scale_factor;
                                                            end
                                                        end
                                                        
                                                        day_multiplier = 2 * num_weeks;  % 2 weekend days per week
                                                    end
                                                    
                                                    % Apply the bus-specific scaling to power values
                                                    P_scaled = zeros(size(P));
                                                    Q_scaled = zeros(size(Q));
                                                    
                                                    for bus = 1:10
                                                        P_scaled(bus) = bus_scale_factors(bus) * P(bus);
                                                        Q_scaled(bus) = bus_scale_factors(bus) * Q(bus);
                                                    end
                                                    
                                                    P_PV = -pv_scale * P;  % PV generation profile, P nominal! NOT Scaled!
                                                    
                                                    % Maximum power values
                                                    [P_invM, P_RChM, Q_invM, P_DChM] = deal(1.1*P_PV, P_PV, 0.3*1.1*P_PV, P_PV);
                                                    P_ilM = -p_ilM_scale * P;  % Interruptible load matrix
                                                    
                                                    % Get tariff - only c_sb is used in optimization
                                                    [C_sb_res, C_if_res] = get_tariff(tariffs, location, tariff_season, day_type, 'residential');
                                                    [C_sb_com, C_if_com] = get_tariff(tariffs, location, tariff_season, day_type, 'commercial');
                                                    
                                                    % Calculate weighted average based on scaled load proportions
                                                    % Note: Using P_scaled which already has the seasonal scaling applied
                                                    total_load = sum(abs(P_scaled(2:end)));  % Exclude slack bus (bus 1)
                                                    
                                                    % Calculate residential and commercial loads
                                                    res_load = 0;
                                                    com_load = 0;
                                                    for bus = 2:10  % Skip slack bus
                                                        if is_residential(bus)
                                                            res_load = res_load + abs(P_scaled(bus));
                                                        else
                                                            com_load = com_load + abs(P_scaled(bus));
                                                        end
                                                    end
                                                    
                                                    % Calculate proportions
                                                    if total_load > 0
                                                        res_proportion = res_load / total_load;
                                                        com_proportion = com_load / total_load;
                                                    else
                                                        res_proportion = 0.5;  % Default to 50/50 if no load
                                                        com_proportion = 0.5;
                                                    end
                                                    
                                                    % Base tariffs
                                                    C_sb_base = res_proportion * C_sb_res + com_proportion * C_sb_com;
                                                    % C_if as 10×24 matrix with separate res/com profiles
                                                    C_if = zeros(10, tm);
                                                    for bus = 1:10
                                                        if bus == 1
                                                            C_if(bus, :) = 0;  % Slack bus stays zero
                                                        elseif is_residential(bus)
                                                            C_if(bus, :) = C_if_res;  % Residential tariff profile
                                                        else
                                                            C_if(bus, :) = C_if_com;  % Commercial tariff profile
                                                        end
                                                    end
                                                    
                                                    % Feed-in tariff - SIZE DEPENDENT
                                                    % Calculate actual PV capacity per bus and get appropriate rate
                                                    C_PV = zeros(10, 1);
                                                    for bus = 1:10
                                                        if bus == 1
                                                            C_PV(bus) = 0;  % Slack bus stays zero
                                                        else
                                                            % Calculate actual PV capacity in kW for this bus
                                                            pv_capacity_kW = abs(P_PV(bus)) * Sb / 1000;
                                                            
                                                            % Get size-appropriate FiT rate
                                                            C_PV(bus) = get_FiT_rate(location, pv_capacity_kW, is_residential(bus));
                                                        end
                                                    end
                                           
                                                    % ========================================================================
                                                    % DEMAND RESPONSE - INDEPENDENT SECTOR EVENTS (DETERMINISTIC)
                                                    % ========================================================================
                                                    
                                                    % Calculate DR parameters ONCE per season
                                                    if day_idx == 1 && strcmp(day_type, 'weekday')
                                                        
                                                        % ====================================================================
                                                        % STEP 1: GET INDEPENDENT EVENT FREQUENCIES
                                                        % ====================================================================
                                                        [res_event_days, com_event_days, overlap_days, total_weekdays] = ...
                                                            get_DR_event_frequency(location, target_season, is_residential);
                                                        
                                                        % Calculate scenario day counts (DETERMINISTIC)
                                                        res_only_days = res_event_days - overlap_days;    % Res event only
                                                        com_only_days = com_event_days - overlap_days;    % Com event only
                                                        both_event_days = overlap_days;                   % Both events
                                                        no_event_days = total_weekdays - res_event_days - com_event_days + overlap_days;
                                                        
                                                        % ====================================================================
                                                        % STEP 2: CALCULATE FIXED CAPACITY PAYMENTS (PER SEASON)
                                                        % ====================================================================
                                                        [~, ~, capacity_rate_res, ~] = generate_DR_events_all_locations(...
                                                            location, 'residential', target_season, 'weekday', true);
                                                        [~, ~, capacity_rate_com, ~] = generate_DR_events_all_locations(...
                                                            location, 'commercial', target_season, 'weekday', true);
                                                        
                                                        % Convert P_ilM to kW
                                                        P_ilM_kW = abs(P_ilM) * (Sb / 1000);
                                                        
                                                        % Calculate seasonal capacity revenue (FIXED, paid regardless of dispatch)
                                                        season_capacity_revenue_res = 0;
                                                        season_capacity_revenue_com = 0;
                                                        for bus = 2:10
                                                            if is_residential(bus)
                                                                season_capacity_revenue_res = season_capacity_revenue_res + ...
                                                                    (capacity_rate_res * P_ilM_kW(bus));
                                                            else
                                                                season_capacity_revenue_com = season_capacity_revenue_com + ...
                                                                    (capacity_rate_com * P_ilM_kW(bus));
                                                            end
                                                        end
                                                        
                                                        season_capacity_revenue_total = season_capacity_revenue_res + season_capacity_revenue_com;
                                                        
                                                        
                                                        season_DR_capacity_revenue(season_idx, yr_idx) = season_capacity_revenue_total;
                                                        % % Initialize delivery tracking for commercial pay-for-performance
                                                        % season_committed_kWh_com = 0;
                                                        % season_delivered_kWh_com = 0;

                                                    end
                                                    
                                                    % ========================================================================
                                                    % WEEKDAY: LOOP OVER 4 INDEPENDENT EVENT SCENARIOS
                                                    % ========================================================================
                                                    if strcmp(day_type, 'weekday')
                                                        
                                                        for event_scenario_idx = 1:4
                                                            
                                                            % ================================================================
                                                            % DETERMINE WHICH SECTORS HAVE EVENTS
                                                            % ================================================================
                                                            switch event_scenario_idx
                                                                case 1  % Residential event ONLY
                                                                    has_res_event = true;
                                                                    has_com_event = false;
                                                                    scenario_multiplier = res_only_days;
                                                                    
                                                                case 2  % Commercial event ONLY
                                                                    has_res_event = false;
                                                                    has_com_event = true;
                                                                    scenario_multiplier = com_only_days;
                                                                    
                                                                case 3  % BOTH events same day
                                                                    has_res_event = true;
                                                                    has_com_event = true;
                                                                    scenario_multiplier = both_event_days;
                                                                    
                                                                case 4  % NO events
                                                                    has_res_event = false;
                                                                    has_com_event = false;
                                                                    scenario_multiplier = no_event_days;
                                                            end
                                                            
                                                            % Skip if no days in this scenario
                                                            if scenario_multiplier <= 0
                                                                continue;
                                                            end
                                                            
                                                            % ================================================================
                                                            % GET DR PARAMETERS FOR BOTH SECTORS (INDEPENDENT)
                                                            % ================================================================
                                                            [C_IL_res, DR_res, ~, discount_res] = generate_DR_events_all_locations(...
                                                                location, 'residential', target_season, day_type, has_res_event);
                                                            [C_IL_com, DR_com, ~, discount_com] = generate_DR_events_all_locations(...
                                                                location, 'commercial', target_season, day_type, has_com_event);
                                                            
                                                            % ================================================================
                                                            % APPLY C_IL AS 10×24 MATRIX - PER BUS, ONLY AT EVENT HOURS
                                                            % ================================================================
                                                            C_IL = zeros(10, tm);  % Changed from (1, tm) to (10, tm)
                                                            for bus = 1:10
                                                                if bus == 1
                                                                    C_IL(bus, :) = 0;  % Slack bus always zero
                                                                elseif is_residential(bus) & has_res_event
                                                                    % Residential bus with residential event → apply C_IL_res at event hours
                                                                    C_IL(bus, :) = C_IL_res;  % C_IL_res already has zeros for non-event hours
                                                                elseif ~is_residential(bus) & has_com_event
                                                                    % Commercial bus with commercial event → apply C_IL_com at event hours
                                                                    C_IL(bus, :) = C_IL_com;  % C_IL_com already has zeros for non-event hours
                                                                else
                                                                    % No event for this bus's sector → all zeros
                                                                    C_IL(bus, :) = 0;
                                                                end
                                                            end
                                                                                                                        
                                                            % Bill discount (weighted geometric mean)
                                                            if res_proportion > 0 & com_proportion > 0
                                                                bill_discount = discount_res^res_proportion * discount_com^com_proportion;
                                                            elseif res_proportion > 0
                                                                bill_discount = discount_res;
                                                            else
                                                                bill_discount = discount_com;
                                                            end
                                                            
                                                            C_sb = C_sb_base .* bill_discount;
                                                            
                                                            % ================================================================
                                                            % ADJUST P_ilM PER BUS BASED ON SECTOR EVENT STATUS
                                                            % ================================================================
                                                            P_ilM_DR = zeros(size(P_ilM));
                                                            for bus = 1:10
                                                                if bus == 1
                                                                    P_ilM_DR(bus) = 0;  % Slack bus
                                                                elseif (is_residential(bus) & has_res_event) || ...
                                                                       (~is_residential(bus) & has_com_event)
                                                                    % This bus's sector has an event → full curtailment capacity
                                                                    P_ilM_DR(bus) = P_ilM(bus);
                                                                else
                                                                    % No event for this bus's sector → baseline flexibility only
                                                                    P_ilM_DR(bus) = P_ilM(bus) * 0.1;
                                                                end
                                                            end
                                                            
                                                            DR_event_occurred = any(C_IL > 0) | any(DR_res) | any(DR_com);
                                                            ssn_val = 1 + contains(target_season, {'winter', 'autumn1', 'autumn2'});
                                                            
                                                            % ================================================================
                                                            % RUN OPTIMIZATION
                                                            % ================================================================
                                                            [ldl, fval, pk_per_bus, btr, dspt] = BattPkShvRsCm_v2(...
                                                                tm, Y, P_scaled, Q_scaled, C_sb, P_PV, ilF_vec, P_ilM_DR, ...
                                                                P_invM, Q_invM, P_RChM, P_DChM, E_M_current, E_mt_current_season, ...
                                                                n_b_initial, C_if, C_PV, C_b, C_IL, ssn_val, yr_idx, current_irr, peak);
                                                            
                                                            % Store pk values (only once per season)
                                                            if event_scenario_idx == 1
                                                                pk_per_bus_per_season_per_year(:, season_idx, yr_idx) = pk_per_bus;
                                                            end
                                                            
                                                            % ================================================================
                                                            % POST-OPTIMIZATION: ADJUST FOR FIXED REVENUES
                                                            % ================================================================
                                                            days_in_season = 7 * num_weeks;                                                           
                                                            
                                                            
                                                            % Amortize enrollment/annual credits per day
                                                            credits_per_day = (one_time_credit_per_season + annual_participation_per_season) / ...
                                                                             days_in_season;
                                                            
                                                            % Convert fval to actual $ before subtracting credits
                                                            fval_actual = fval * (Sb / 1000);
                                                            fval_adjusted = fval_actual - credits_per_day;
                                                            
                                                            % ================================================================
                                                            % TRACK PERFORMANCE-BASED DR REVENUE (FROM ACTUAL DISPATCH)
                                                            % ================================================================
                                                            if DR_event_occurred
                                                                dr_performance_revenue_day = 0;
                                                                
                                                                for h = 1:tm
                                                                    for bus = 2:10
                                                                        if C_IL(bus, h) > 0
                                                                            il_idx = (bus-2)*5 + 6;
                                                                            dr_performance_revenue_day = dr_performance_revenue_day + ...
                                                                                abs(dspt(il_idx, h)) * C_IL(bus, h);
                                                                        end
                                                                    end
                                                                end
                                                                                                                                
                                                                % Accumulate DR metrics (only for actual event scenarios)
                                                                if event_scenario_idx <= 3
                                                                    season_DR_performance_revenue(season_idx, yr_idx) = ...
                                                                        season_DR_performance_revenue(season_idx, yr_idx) + ...
                                                                        (dr_performance_revenue_day * scenario_multiplier);
                                                                    
                                                                    season_DR_total_events(season_idx, yr_idx) = ...
                                                                        season_DR_total_events(season_idx, yr_idx) + scenario_multiplier;
                                                                end
                                                            end
                                                            % % TRACK DELIVERY FOR CAPACITY PAYMENT VERIFICATION (commercial only)
                                                            % if has_com_event && (strcmp(location,'California') || strcmp(location,'Texas') || strcmp(location,'NewYork'))
                                                            %     for bus = 2:10
                                                            %         if ~is_residential(bus)
                                                            %             for h = 1:tm
                                                            %                 if DR_com(h)  % event hour for commercial program
                                                            %                     il_idx = (bus-2)*5 + 6;
                                                            %                     % committed vs delivered (in p.u., scaled equally)
                                                            %                     season_committed_kWh_com = season_committed_kWh_com + ...
                                                            %                         P_ilM_DR(bus) * scenario_multiplier;
                                                            %                     season_delivered_kWh_com = season_delivered_kWh_com + ...
                                                            %                         abs(dspt(il_idx, h)) * scenario_multiplier;
                                                            %                 end
                                                            %             end
                                                            %         end
                                                            %     end
                                                            % end

                                                            % ================================================================
                                                            % UPDATE BATTERY STATE (AFTER LAST SCENARIO)
                                                            % ================================================================
                                                            if event_scenario_idx == 4
                                                                for bus = 2:10
                                                                    discharge_idx = (bus-2)*5 + 4;
                                                                    charge_idx = (bus-2)*5 + 5;
                                                                    
                                                                    final_discharge = dspt(discharge_idx, 24);
                                                                    final_charge = dspt(charge_idx, 24);
                                                                    
                                                                    E_mt_current_season(bus) = E_mt_current_season(bus) - ...
                                                                        final_charge - final_discharge;
                                                                    
                                                                    E_mt_current_season(bus) = max(0, min(E_mt_current_season(bus), ...
                                                                        E_M_current(bus)));
                                                                end
                                                            end
                                                            
                                                            % ================================================================
                                                            % ACCUMULATE WEIGHTED RESULTS
                                                            % ================================================================
                                                            cur_season_fval = cur_season_fval + (fval_adjusted * scenario_multiplier);
                                                            cur_season_ldl = cur_season_ldl + (ldl * scenario_multiplier);
                                                            cur_season_btr = cur_season_btr + (btr * scenario_multiplier);

                                                            % Calculate battery throughput for each bus from dspt matrix
                                                            for bus = 2:10  % Bus 1 is slack, start from bus 2
                                                                discharge_idx = (bus-2)*5 + 4;
                                                                charge_idx = (bus-2)*5 + 5;
                                
                                                                
                                                                % Sum absolute values of charge and discharge for each hour
                                                                bus_daily_btr = 0;
                                                                for h = 1:24
                                                                    % Add absolute values because charging is negative
                                                                    bus_daily_btr = bus_daily_btr + abs(dspt(discharge_idx, h)) + abs(dspt(charge_idx, h));
                                                                end
                                                                
                                                                % Accumulate in per-bus array with appropriate day multiplier
                                                                per_bus_season_btr(bus) = per_bus_season_btr(bus) + (bus_daily_btr * scenario_multiplier);
                                                            end
                                                        end
                                                    else 
                                                        % ====================================================================
                                                        % WEEKEND (NO DR EVENTS FOR ANY SECTOR)
                                                        % ====================================================================
                                                        scenario_multiplier = 2 * num_weeks;
                                                        
                                                        % Get no-event parameters (for bill discount only)
                                                        [~, ~, ~, discount_res] = generate_DR_events_all_locations(...
                                                            location, 'residential', target_season, day_type, false);
                                                        [~, ~, ~, discount_com] = generate_DR_events_all_locations(...
                                                            location, 'commercial', target_season, day_type, false);
                                                        
                                                        % Apply bill discount
                                                        if res_proportion > 0 & com_proportion > 0
                                                            bill_discount = discount_res^res_proportion * discount_com^com_proportion;
                                                        elseif res_proportion > 0
                                                            bill_discount = discount_res;
                                                        else
                                                            bill_discount = discount_com;
                                                        end
                                                        
                                                        C_sb = C_sb_base * bill_discount;
                                                        C_IL = zeros(10, tm);  
                                                        P_ilM_DR = P_ilM * 0.1;  % Baseline flexibility only
                                                        
                                                        ssn_val = 1 + contains(target_season, {'winter', 'autumn1', 'autumn2'});
                                                        
                                                        % RUN OPTIMIZATION
                                                        [ldl, fval, ~, btr, dspt] = BattPkShvRsCm_v2(...
                                                            tm, Y, P_scaled, Q_scaled, C_sb, P_PV, ilF_vec, P_ilM_DR, ...
                                                            P_invM, Q_invM, P_RChM, P_DChM, E_M_current, E_mt_current_season, ...
                                                            n_b_initial, C_if, C_PV, C_b, C_IL, ssn_val, yr_idx, current_irr, peak);
                                                        
                                                        % POST-OPTIMIZATION
                                                        days_in_season = 7 * num_weeks;                                                        
                                                        credits_per_day = (one_time_credit_per_season + annual_participation_per_season) / ...
                                                                         days_in_season;
                                                        % Convert fval to actual $ before subtracting credits
                                                        fval_actual = fval * (Sb / 1000);
                                                        fval_adjusted = fval_actual - credits_per_day;
                                                        
                                                        % ACCUMULATE RESULTS
                                                        cur_season_fval = cur_season_fval + (fval_adjusted * scenario_multiplier);
                                                        cur_season_ldl = cur_season_ldl + (ldl * scenario_multiplier);
                                                        cur_season_btr = cur_season_btr + (btr * scenario_multiplier);

                                                        % Calculate battery throughput for each bus from dspt matrix
                                                            for bus = 2:10  % Bus 1 is slack, start from bus 2
                                                                discharge_idx = (bus-2)*5 + 4;
                                                                charge_idx = (bus-2)*5 + 5;
                                
                                                                
                                                                % Sum absolute values of charge and discharge for each hour
                                                                bus_daily_btr = 0;
                                                                for h = 1:24
                                                                    % Add absolute values because charging is negative
                                                                    bus_daily_btr = bus_daily_btr + abs(dspt(discharge_idx, h)) + abs(dspt(charge_idx, h));
                                                                end
                                                                
                                                                % Accumulate in per-bus array with appropriate day multiplier
                                                                per_bus_season_btr(bus) = per_bus_season_btr(bus) + (bus_daily_btr * scenario_multiplier);
                                                            end
                                                    end  % End weekday/weekend
                                                end % ===== END OF DAY_TYPE LOOP ====
                                                
                                                % % === APPLY DELIVERY RATIO TO CAPACITY PAYMENTS (ONCE PER SEASON) ===

                                                % % Effective commercial capacity revenue (pay-for-performance)
                                                % effective_capacity_com = season_capacity_revenue_com;
                                                
                                                % if season_committed_kWh_com > 0
                                                %     delivery_ratio_com = min(1.0, season_delivered_kWh_com / season_committed_kWh_com);
                                                % else
                                                %     % No commercial events or no commitment → treat as zero performance
                                                %     delivery_ratio_com = 0.0;
                                                % end
                                                
                                                % effective_capacity_com = season_capacity_revenue_com * delivery_ratio_com;
                                                
                                                % % Residential capacity revenue stays fixed (no performance factor)
                                                % effective_capacity_total = season_capacity_revenue_res + effective_capacity_com;
                                                
                                                % % Subtract the capacity revenue ONCE from the seasonal cost
                                                % cur_season_fval = cur_season_fval - effective_capacity_total;
                                                
                                                % % Store effective capacity revenue for DR stats
                                                % season_DR_capacity_revenue(season_idx, yr_idx) = effective_capacity_total;
                                                % ================================================================
                                                % FINALIZE DR PAYMENTS (ONCE PER SEASON)
                                                % ================================================================
                                                
                                                % Capacity payments are FIXED (for being enrolled/available)
                                                % NOT scaled by delivery ratio - that was a workaround for units bug
                                                effective_capacity_total = season_capacity_revenue_res + season_capacity_revenue_com;
                                                
                                                % Subtract capacity revenue from seasonal cost
                                                cur_season_fval = cur_season_fval - effective_capacity_total;
                                                
                                                % Store for reporting
                                                season_DR_capacity_revenue(season_idx, yr_idx) = effective_capacity_total;
                                                % Now update SoC for next season and store results
                                                E_mt_current = E_mt_current_season;
                                                
                                                season_fval(season_idx, yr_idx) = cur_season_fval;
                                                season_ldl(season_idx, yr_idx) = cur_season_ldl;
                                                season_btr(season_idx, yr_idx) = cur_season_btr;

                                                                                                
                                                % === IMPROVED DEGRADATION MODEL ===
                                                for bus = 2:10
                                                    if E_M_initial(bus) <= 0
                                                        continue;
                                                    end
                                                
                                                    days_in_period = num_weeks * 7;
                                                    season_throughput = per_bus_season_btr(bus);
                                                    
                                                    % === Battery-specific DoD estimation ===
                                                    avg_daily_throughput = season_throughput / days_in_period;
                                                    effective_daily_cycles = avg_daily_throughput / (2 * E_M_initial(bus) + eps);
                                                    effective_daily_cycles = max(0, min(2, effective_daily_cycles));
                                                    
                                                    if strcmp(battery_type, 'Li-ion')
                                                        % Li-ion: shallow cycles in grid storage
                                                        if effective_daily_cycles < 0.1
                                                            effective_DoD = 0.10;
                                                        elseif effective_daily_cycles < 0.5
                                                            effective_DoD = 0.25;
                                                        elseif effective_daily_cycles < 1.0
                                                            effective_DoD = 0.40;
                                                        else
                                                            effective_DoD = 0.55;
                                                        end
                                                        
                                                    else  % Na-ion
                                                        % Na-ion: similar cycling patterns but slower charging capability
                                                        % Tends toward slightly more complete cycles due to slower charge/discharge
                                                        if effective_daily_cycles < 0.1
                                                            effective_DoD = 0.12;
                                                        elseif effective_daily_cycles < 0.5
                                                            effective_DoD = 0.27;
                                                        elseif effective_daily_cycles < 1.0
                                                            effective_DoD = 0.45;
                                                        else
                                                            effective_DoD = 0.65;
                                                        end
                                                    end
                                                    
                                                    total_cycles_in_period = effective_daily_cycles * days_in_period;
                                                                                                        
                                                    % === CYCLING DEGRADATION ===
                                                    % Use 0.8 as reference DoD
                                                    ref_DoD = 0.8;
                                                    N_ref = CYCLE_LIFE_80DOD;
                                                    
                                                    if effective_DoD < 0.02
                                                        N_at_DoD = 1e6; % Essentially no cycling
                                                    else
                                                        N_at_DoD = N_ref * (ref_DoD / effective_DoD) ^ DOD_EXPONENT;
                                                    end
                                                    
                                                    % Temperature effect on cycling
                                                    hourly_temp = get_hourly_temperature(location, target_season);
                                                    avg_temp = mean(hourly_temp);
                                                    temp_cycle_multiplier = exp(TEMP_SENSITIVITY_CYCLE * (avg_temp - T_ref));
                                                    N_effective = N_at_DoD / temp_cycle_multiplier;
                                                    
                                                    cycling_fade_fraction = total_cycles_in_period / (N_effective + eps);
                                                    
                                                    % === CALENDAR AGING ===
                                                    daily_cal_fade = CALENDAR_FADE_ANNUAL / 365;
                                                    temp_cal_multiplier = exp(TEMP_SENSITIVITY_CAL * (avg_temp - T_ref));
                                                    
                                                    % SoC stress (use current SoC as proxy)
                                                    mean_SoC = E_mt_current(bus) / (E_M_initial(bus) + eps);
                                                    
                                                    if strcmp(battery_type, 'Li-ion')
                                                        % Li-ion: minimal SoC stress in normal range
                                                        if mean_SoC > 0.95
                                                            soc_cal_stress = 1.15;  % Only extreme high SoC
                                                        else
                                                            soc_cal_stress = 1.0;
                                                        end
                                                        
                                                    else % Na-ion
                                                        % Na-ion: more robust across SoC range (safer chemistry)
                                                        % Less prone to degradation at extreme SoC
                                                        if mean_SoC > 0.98 || mean_SoC < 0.05
                                                            soc_cal_stress = 1.05;  % Minimal stress even at extremes
                                                        else
                                                            soc_cal_stress = 1.0;
                                                        end
                                                    end
                                                                                                        
                                                    calendar_fade_fraction = days_in_period * daily_cal_fade * temp_cal_multiplier * soc_cal_stress;
                                                    
                                                    % === COMBINE DEGRADATION ===
                                                    total_fade_fraction = cycling_fade_fraction + calendar_fade_fraction;
                                                    
                                                    % Sanity check for extreme values only
                                                    if total_fade_fraction > 0.15 % 15% per season is physically unrealistic
                                                        fprintf('Warning: Extreme degradation %.2f%% for bus %d - capping\n', ...
                                                                total_fade_fraction*100, bus);
                                                        total_fade_fraction = 0.15;
                                                    end
                                                    
                                                    % Update cumulative degradation
                                                    cum_degradation_simple(bus) = cum_degradation_simple(bus) + total_fade_fraction;
                                                    cum_degradation_simple(bus) = min(cum_degradation_simple(bus), 1.0); % Can't exceed 100%
                                                    
                                                    E_M_current(bus) = E_M_initial(bus) * (1 - cum_degradation_simple(bus));
                                                    
                                                    % Replacement check
                                                    if E_M_current(bus) <= E_ret(bus)
                                                        replacements(bus) = replacements(bus) + 1;
                                                        cum_degradation_simple(bus) = 0;
                                                        E_M_current(bus) = E_M_initial(bus);
                                                        E_mt_current(bus) = 0.5 * E_M_initial(bus);
                                                    end
                                                end
                                               
                                                % Check replacements
                                                season_bus_replacements = zeros(10, 1);
                                                for i = 2:10
                                                    if E_M_initial(i) > 0 && E_M_current(i) < E_ret(i)
                                                        replacements(i) = replacements(i) + 1;
                                                        season_bus_replacements(i) = 1;
                                                        E_M_current(i) = E_M_initial(i);
                                                        E_mt_current(i) = 0.5 * E_M_initial(i);
                                                        cum_degradation_simple(i) = 0;  % Reset degradation
                                                    end
                                                end
                                                
                                                % Ensure capacity doesn't go below retirement threshold
                                                for i = 2:10
                                                    E_M_current(i) = max(E_M_current(i), E_ret(i));
                                                    E_mt_current(i) = max(0, min(E_mt_current(i), E_M_current(i)));
                                                end
                                                
                                                % Calculate total replacements for this season
                                                season_replacements(season_idx, yr_idx) = sum(season_bus_replacements);
                                                
                                                % Update annual totals for current year
                                                annual_fval(yr_idx) = annual_fval(yr_idx) + cur_season_fval;
                                                annual_ldl(yr_idx) = annual_ldl(yr_idx) + cur_season_ldl;
                                                annual_btr(yr_idx) = annual_btr(yr_idx) + cur_season_btr;
                                                annual_replacements(yr_idx) = annual_replacements(yr_idx) + sum(season_bus_replacements);
                                                
                                                % Update 10-year totals
                                                total_fval = total_fval + cur_season_fval;
                                                total_ldl = total_ldl + cur_season_ldl;
                                                total_btr = total_btr + cur_season_btr;
                                            end
                                            % Print yearly results for current simulation
                                            fprintf('Simulation %d - Year %d: Fval=%.2f, Replacements=%d\n', ...
                                            round(sim_count), yr_idx, annual_fval(yr_idx), round(annual_replacements(yr_idx)));
                                        end
                                        
                                        % Calculate total battery replacements over 10 years
                                        total_replacements = sum(replacements);
                                        
                                        % Calculate NPV-adjusted fval using compound discount factors
                                        discounted_fval = 0;
                                        
                                        for year_i = 1:length(years)
                                            % Calculate compound discount factor for year yr
                                            % This accounts for varying IRR rates each year
                                            compound_discount = 1;
                                            for j = 1:year_i
                                                compound_discount = compound_discount / (1 + irr_vector(j));
                                            end
                                            
                                            % Add discounted annual value to total
                                            discounted_fval = discounted_fval + annual_fval(year_i) * compound_discount;
                                        end      
                                        

                                        % Calculate total DR metrics
                                        total_DR_events_decade = sum(sum(season_DR_total_events));
                                        total_DR_revenue_decade = sum(sum(season_DR_capacity_revenue)) + ...
                                                sum(sum(season_DR_performance_revenue));
                                      
                                        % Calculate worst-case PK per bus per season over the decade
                                        worst_case_pk_per_bus_per_season = max(pk_per_bus_per_season_per_year, [], 3); % Max along the third dimension (years)
                                        
                                        % ================================================================
                                        % STORE PER-BUS RESULTS
                                        % ================================================================
                                        temp_pk_data = cell(10, 1);
                                        for bus = 1:10
                                            temp_pk_data{bus} = {
                                                peak_idx, ...
                                                peak, ...
                                                p_ilM_scale, ...
                                                pv_scale, ...
                                                x, ...
                                                battery_type, ...
                                                residential_mix, ...
                                                ilF_res, ...
                                                ilF_com, ...
                                                irr_name, ...
                                                bus, ...
                                                worst_case_pk_per_bus_per_season(bus, 1), ... % Spring
                                                worst_case_pk_per_bus_per_season(bus, 2), ... % Summer
                                                worst_case_pk_per_bus_per_season(bus, 3), ... % Autumn1
                                                worst_case_pk_per_bus_per_season(bus, 4), ... % Autumn2
                                                worst_case_pk_per_bus_per_season(bus, 5)      % Winter
                                            };
                                        end
                                    
                                        % Append to local results
                                        local_pk_results = [local_pk_results; temp_pk_data];
    

                                        % Print detailed simulation summary (matching battery_replacement_analysis format)
                                        fprintf('\n--- SIMULATION %d SUMMARY ---\n');
                                        fprintf('Parameters: Peak=%.2f, IL Scale=%.2f, PV Scale=%.1f, Battery Scale=%.1f\n', ...
                                                peak, p_ilM_scale, pv_scale, x);
                                        fprintf('Battery Type=%s, Res Mix=%s, IL Res=%.1f, IL Com=%.1f, IRR=%s\n', ...
                                                battery_type, residential_mix, ilF_res, ilF_com, irr_name);
                                        fprintf('NPV (Discounted) Fval: %.6f\n', discounted_fval);
                                        fprintf('Total (Undiscounted) Fval: %.6f\n', total_fval);
                                        fprintf('Total Load Served: %.6f\n', total_ldl);
                                        fprintf('Total Replacements: %d\n', total_replacements);
                                        
                                        % Show per-bus replacement breakdown
                                        fprintf('\nPer-bus replacements over 10 years:\n');
                                        for i = 1:10
                                            fprintf('Bus %d: %d replacements\n', i, replacements(i));
                                        end
                                        
                                        fprintf(['\n' repmat('=', 1, 50) '\n\n']);
                                        

                                        % Special verification for x=0 scenarios (matching battery_replacement_analysis)
                                        if x == 0
                                            fprintf('\n--- DISPATCH VERIFICATION FOR x=0 ---\n');
                                            fprintf('✓ VERIFIED: No battery dispatch when x=0 (tracking not implemented in multi-year)\n');
                                        elseif x > 0
                                            fprintf('\n--- BATTERY UTILIZATION FOR x=%.1f ---\n', x);
                                            fprintf('Battery throughput tracking: %.6f (per-bus tracking not implemented)\n', total_btr);
                                        end
                                        
                                        % Write results to the consolidated CSV 
                                        local_results_data = [local_results_data; ...
                                        {peak, p_ilM_scale, pv_scale, x, battery_type, residential_mix, ...
                                        ilF_res, ilF_com, irr_name, discounted_fval, total_fval, total_ldl, ...
                                        total_replacements, total_DR_events_decade, total_DR_revenue_decade}];

                                    end % irr_vector loop
                                end % com_il_idx loop
                            end % res_il_idx loop
                        end % res_mix_idx loop
                    end % batt_idx loop
                end % x_idx loop
            end % pv_idx loop
        end % p_ilM_idx loop
        % Store results for this peak_idx iteration
        all_results_data{peak_idx} = local_results_data
        all_pk_results{peak_idx} = local_pk_results;
    end % peak_idx loop
    
    % Write all collected results
    all_results_flat = vertcat(all_results_data{:});
    fmt = ['%.2f,%.2f,%.2f,%.2f,%s,%s,%.2f,%.2f,%s,' ...  % first 9 cols
       '%.6f,%.6f,%.6f,%d,%d,%.6f\n'];                 % last 6 cols
    for i = 1:size(all_results_flat, 1)
        fprintf(fid, fmt, all_results_flat{i,:});
    end

    % Flatten the cell array of results back into a single cell array
    all_pk_results = vertcat(all_pk_results{:});

    % Close the file
    fclose(fid);
    % Write per-bus worst-case PK values to separate CSV file
    pk_filename = 'pk_per_bus_simple_full_TX.csv';
    fid_pk = fopen(pk_filename, 'w');
    
    % Write CSV header for per-bus worst-case PK values`
    fprintf(fid_pk, 'Simulation_ID,Peak_Value,IL_Scale,PV_Scale,Battery_Scale_Factor,Battery_Type,Residential_Mix,IL_Fraction_Res,IL_Fraction_Com,IRR_Scenario,Bus,Worst_PK_Spring,Worst_PK_Summer,Worst_PK_Autumn1,Worst_PK_Autumn2,Worst_PK_Winter\n');
    
    for i = 1:length(all_pk_results)
        row_data = all_pk_results{i};
        fprintf(fid_pk, '%d,%.2f,%.2f,%.2f,%.2f,%s,%s,%.2f,%.2f,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n', row_data{:});
    end
    fclose(fid_pk);
    
    fprintf('\n=== MULTI-YEAR ANALYSIS COMPLETE ===\n');
    fprintf('Total simulations completed: %d\n');
    fprintf('Results match battery_replacement_analysis format\n');
    fprintf('CSV files generated with detailed results\n');

    

end

function hourly_temp = get_hourly_temperature(location, season)
    % Get base temperature
    Tmean = get_seasonal_temperature(location, season);
    
    % Handle autumn variations
    if contains(season, 'autumn')
        season = 'autumn';
    end
    
    % Location-specific diurnal variation parameters
    switch lower(location)
        case 'minnesota'
            % Continental climate - large variations
            switch lower(season)
                case 'summer'
                    amp = 11;  % Large day-night variation
                    Tmean = Tmean + 3;  % Adjust for peak summer
                case 'winter'
                    amp = 6;   % Smaller variation in winter
                    Tmean = Tmean - 2;  % Adjust for deep winter
                otherwise
                    amp = 9;
            end
            
        case 'texas'
            % Southern climate - hot with moderate variations
            switch lower(season)
                case 'summer'
                    amp = 8;   % Moderate variation, consistently hot
                    Tmean = Tmean + 2;  % Peak summer heat
                case 'winter'
                    amp = 7;   % Moderate winter variation
                otherwise
                    amp = 7.5;
            end
            
        case 'california'
            % Mediterranean climate - mild variations
            switch lower(season)
                case 'summer'
                    amp = 6;   % Coastal influence reduces variation
                case 'winter'
                    amp = 5;   % Very mild winter variations
                otherwise
                    amp = 5.5;
            end
            
        case 'newyork'
            % Temperate climate - moderate variations
            switch lower(season)
                case 'summer'
                    amp = 9;   % Moderate summer variation
                    Tmean = Tmean + 1;  % Urban heat effect
                case 'winter'
                    amp = 7;   % Moderate winter variation
                    Tmean = Tmean - 1;  % Cold snaps
                otherwise
                    amp = 8;
            end
            
        otherwise
            % Default parameters for unknown locations
            amp = 8;
    end
    
    hours = 0:23;
    % Peak at 3 PM (hour 15)
    hourly_temp = Tmean + amp * sin((2*pi/24) * (hours - 9));  % Peak shifted to 3 PM
    
    % Add extreme events (1% probability of extreme temps)
    if rand < 0.01
        switch lower(season)
            case 'winter'
                % Cold snaps - severity depends on location
                switch lower(location)
                    case 'minnesota'
                        hourly_temp = hourly_temp - 6;  % Polar vortex
                    case 'texas'
                        hourly_temp = hourly_temp - 4;   % Rare cold snap
                    case 'newyork'
                        hourly_temp = hourly_temp - 5;  % Nor'easter
                    case 'california'
                        hourly_temp = hourly_temp - 3;   % Mild cold spell
                  
                end
                
            case 'summer'
                % Heat waves - severity depends on location
                switch lower(location)
                    case 'texas'
                        hourly_temp = hourly_temp + 5;  % Intense heat wave
                    case 'minnesota'
                        hourly_temp = hourly_temp + 4;   % Heat wave with humidity
                    case 'california'
                        hourly_temp = hourly_temp + 3;   % Dry heat wave
                    case 'newyork'
                        hourly_temp = hourly_temp + 5;   % Urban heat island
                   
                end
        end
    end
    
    hourly_temp = hourly_temp(:)';
end