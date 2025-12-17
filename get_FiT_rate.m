function C_PV_bus = get_FiT_rate(location, pv_capacity_kW, is_residential)
% GET_FIT_RATE Returns Feed-in Tariff rate based on location, system size, and sector
%
% Based on actual program data:
% - California: FIT for <3MW at $0.08923/kWh (10-20 year contracts)
% - Minnesota: Made in Minnesota program (varies by manufacturer, using averages)
% - Texas: Austin Energy PBI program (tiered by size)
% - New York: Competitive bidding (using historical averages)
%
% Inputs:
%   location        - String: 'California', 'NewYork', 'Texas', 'Minnesota'
%   pv_capacity_kW  - Numeric: PV system capacity in kW
%   is_residential  - Boolean: true for residential/small commercial, false for commercial
%
% Output:
%   C_PV_bus        - Feed-in tariff rate in $/kWh

    switch location
        case 'California'
            % California FIT Program: $0.08923/kWh for systems <3MW
            % Same rate regardless of size (up to 3MW)
            if pv_capacity_kW <= 3000  % 3 MW limit
                C_PV_bus = 0.08923;
            else
                C_PV_bus = 0.05;  % Wholesale rate for larger systems
            end
            
        case 'NewYork'
            % NY-Sun Competitive PV Program
            % Rates based on competitive bidding - using historical averages
            % Strategic location bonus (+15%) not included here
            if pv_capacity_kW < 200
                % Below program minimum - use VDER rates
                if is_residential
                    C_PV_bus = 0.085;  % Small DG VDER
                else
                    C_PV_bus = 0.075;  % Small commercial VDER
                end
            elseif pv_capacity_kW < 500
                % Medium commercial (200-500 kW)
                C_PV_bus = 0.065;  % Competitive bid average
            elseif pv_capacity_kW < 1000
                % Large commercial (500-1000 kW)
                C_PV_bus = 0.055;  % Competitive bid average
            else
                % Very large (>1 MW)
                C_PV_bus = 0.045;  % Lower competitive bids
            end
            
        case 'Texas'
            % Austin Energy PBI Program
            % Tiered by system size
            if is_residential
                % Residential not covered by Austin Energy commercial program
                % Use PEC buyback rate
                C_PV_bus = 0.0826;  % PEC Sustainable Power Credit (March 2025)
            else
                % Commercial tiers per Austin Energy
                if pv_capacity_kW < 200
                    % Small/Medium Commercial <200 kW-AC
                    C_PV_bus = 0.08;
                elseif pv_capacity_kW < 1000
                    % Large Commercial 200-999 kW-AC
                    C_PV_bus = 0.06;
                else
                    % 1+ MW
                    C_PV_bus = 0.04;
                end
            end
            
        case 'Minnesota'
            % Made in Minnesota (MiM) Solar Incentive Program
            % Using average rates across module manufacturers
            % (Silicon Energy, tenKsolar, Heliene, itek Energy)
            if pv_capacity_kW < 10 && is_residential
                % Residential <10 kW DC
                % Average of: $0.30, $0.23, $0.23, $0.27 = $0.2575
                C_PV_bus = 0.2575;
            elseif pv_capacity_kW < 40
                if is_residential
                    % Small residential/non-profit <40 kW
                    % Average of: $0.25, $0.15, $0.15, $0.20 = $0.1875
                    C_PV_bus = 0.1875;
                else
                    % Commercial/For Profit <40 kW
                    % Average of: $0.23, $0.13, $0.13, $0.18 = $0.1675
                    C_PV_bus = 0.1675;
                end
            else
                % Systems >40 kW don't qualify for MiM
                % Fall back to standard net metering / Value of Solar
                if pv_capacity_kW < 500
                    C_PV_bus = 0.08;   % Net metering rate
                else
                    C_PV_bus = 0.065;  % Near-wholesale for large systems
                end
            end
            
        otherwise
            % Default fallback
            warning('Unknown location: %s. Using default rate.', location);
            C_PV_bus = 0.06;
    end
    
    % Ensure non-negative
    C_PV_bus = max(0, C_PV_bus);
end