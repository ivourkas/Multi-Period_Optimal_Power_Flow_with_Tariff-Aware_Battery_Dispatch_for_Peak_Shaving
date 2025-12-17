function [E_M, E_mt, n_b, C_b] = get_battery_characteristics(battery_type, P, x)

    y = 0.5;  % SOC initialization factor
    
    switch battery_type
        case 'Li-ion'
            % Modern LFP grid storage (2024 pricing)
            n_b = 0.92 * ones(size(P));  % 92% round-trip efficiency
            % Cost: $600/kWh, Cycle life: 4500 at 80% DoD
            C_b = 600 / 4500;  
            
        case 'Na-ion'
            % Sodium-ion for grid storage
            n_b = 0.88 * ones(size(P));  % 88% efficiency (slightly lower due to internal resistance)
            % Cost: $400/kWh (reflects 50x cheaper raw materials but lower energy density)
            C_b = 300 / 3000;  
    end
    
    if x == 0
        E_M = zeros(size(P));
        E_mt = zeros(size(P));
    else
        E_M = -x * P;
        E_mt = -y * x * P;
    end
end