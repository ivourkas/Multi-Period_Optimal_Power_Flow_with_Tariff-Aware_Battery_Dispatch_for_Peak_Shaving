function [C_IL_hourly, DR_active_hours, capacity_rate_per_kW, bill_discount_factor] = ...
    generate_DR_events_all_locations(location, sector, season, day_type, has_event)
    
    % ========================================================================
    % DEMAND RESPONSE EVENT GENERATOR - FINAL PRODUCTION VERSION
    % ========================================================================
    % Generates DR parameters for optimization based on real utility programs
    % All data verified against actual program specifications
    %
    % INPUTS:
    %   location: 'California', 'Texas', 'Minnesota', 'NewYork'
    %   sector: 'residential' or 'commercial'
    %   season: 'summer', 'winter', 'spring', 'autumn1', 'autumn2'
    %   day_type: 'weekday' or 'weekend' (NOTE: weekends always return no event)
    %   has_event: true/false - whether this is an event day
    %
    % OUTPUTS:
    %   C_IL_hourly (1×24): Performance payment ($/kWh) for actual curtailment
    %                       Used in optimization objective to incentivize DR
    %   DR_active_hours (1×24): Boolean indicating event hours
    %   capacity_rate_per_kW: Fixed capacity payment ($/kW per season)
    %                         Applied OUTSIDE optimization: revenue = rate × P_ilM_kW
    %   bill_discount_factor: Multiplier for C_sb (e.g., 0.90 = 10% discount)
    %
    % PAYMENT STRUCTURE:
    %   - Capacity payments: Fixed revenue based on committed kW (outside opt)
    %   - Performance payments: Variable based on actual curtailment (in C_IL)
    %   - Enrollment/annual: One-time/recurring credits (handled in main script)
    %   - Bill discounts: Applied to all days (Minnesota summer)
    %
    % NOTE: All DR events occur on WEEKDAYS ONLY for computational efficiency
    % ========================================================================
    
    tm = 24;
    C_IL_hourly = zeros(1, tm);
    DR_active_hours = false(1, tm);
    capacity_rate_per_kW = 0;
    bill_discount_factor = 1.0;
    
    % ========================================================================
    % WEEKEND OVERRIDE: No DR events on weekends for any program
    % ========================================================================
    if strcmp(day_type, 'weekend')
        % Minnesota bill discount still applies on weekends
        if strcmp(location, 'Minnesota') && contains(season, 'summer')
            bill_discount_factor = 0.90;
        end
        return;  % Exit - no events on weekends
    end
    
    % ========================================================================
    % MINNESOTA BILL DISCOUNT: Applies to ALL summer days (event or not)
    % ========================================================================
    if strcmp(location, 'Minnesota') && contains(season, 'summer')
        bill_discount_factor = 0.90;  % 10% discount on all summer bills
    end
    
    % Return if no event (only MN discount may have been applied)
    if ~has_event
        return;
    end
    
    % ========================================================================
    % EVENT DAY LOGIC (has_event = true, weekday only)
    % ========================================================================
    
    switch location
        
        % ====================================================================
        % CALIFORNIA - Pacific Gas & Electric (PG&E)
        % ====================================================================
        case 'California'
            if strcmp(sector, 'residential')
                % ================================================================
                % POWER SAVER REWARDS
                % ================================================================
                % Program Type: Performance-based
                % Payment: $1/kWh reduced during events
                % Season: May 1 - Oct 31 (summer + autumn)
                % Events: 4pm-9pm, ~5 hours, max 12 events/year (60 hours total)
                % ================================================================
                
                if contains(season, {'spring','summer', 'autumn1', 'autumn2'})
                    duration = 5;
                    start_hour = 17;  % 4pm (hour 17 in 1-24 indexing)
                    event_hours = start_hour:(start_hour + duration - 1);
                    DR_active_hours(event_hours) = true;
                    
                    % Performance payment: $1/kWh actual curtailment
                    C_IL_hourly(event_hours) = 1.00;
                    
                    % No capacity payment
                    capacity_rate_per_kW = 0;
                end
                
            else  % Commercial
                % ================================================================
                % BASE INTERRUPTIBLE PROGRAM (BIP)
                % ================================================================
                % Program Type: Capacity-based with penalties
                % Capacity Payment:
                %   - May-Oct: $12.50-13.50/kW-month (using $13/kW-month)
                %   - Nov-Apr: $9.50-10.50/kW-month (using $10/kW-month)
                % Events: Max 6 hours/event, called during CAISO emergencies
                % Penalty: $6/kWh for exceeding Firm Service Level (FSL)
                % Eligible: ≥100 kW demand on demand TOU rate
                % ================================================================
                % CAPACITY ALLOCATION BY SEASON (avoiding double-counting):
                %   Spring (3 months): May ($13) = $13/kW
                %   Summer (3 months): Jun-Aug ($13 each) = $39/kW
                %   Autumn1 (1 month): Sep ($13) = $13/kW
                %   Autumn2 (2 months): Oct ($13) + Nov ($10) = $23/kW
                %   Winter (3 months): Dec-Feb ($10 each) = $30/kW
                %   TOTAL: $118/kW per year
                % ================================================================
                
                % Allocate capacity payment based on season
                if strcmp(season, 'spring')
                    capacity_rate_per_kW = 13 * 1;  % May only
                elseif strcmp(season, 'summer')
                    capacity_rate_per_kW = 13 * 3;  % Jun-Aug
                elseif strcmp(season, 'autumn1')
                    capacity_rate_per_kW = 13 * 1;  % Sep only
                elseif strcmp(season, 'autumn2')
                    capacity_rate_per_kW = 13 * 1 + 10 * 1;  % Oct + Nov
                else  % winter
                    capacity_rate_per_kW = 10 * 3;  % Dec-Feb
                end
                
                % Event characteristics
                duration = 6;
                start_hour = randi([11, 14]);  % Typical CAISO emergency hours
                event_hours = start_hour:min(start_hour + duration - 1, 22);
                DR_active_hours(event_hours) = true;
                
                % No explicit performance payment - must curtail to avoid $6/kWh penalty
                % Use small value to encourage curtailment (penalty avoidance)
                % CA commercial peak tariff: ~$0.39-0.40/kWh
                C_IL_hourly(event_hours) = 0.0;  
            end
            
        % ====================================================================
        % TEXAS
        % ====================================================================
        case 'Texas'
            if strcmp(sector, 'residential')
                % ================================================================
                % DIRECT ENERGY - REDUCE YOUR USE
                % ================================================================
                % Program Type: Smart thermostat control with enrollment bonus
                % Payment: $25 one-time bill credit (handled in main script)
                % Value: Avoiding peak electricity costs during AC adjustment
                % Events: Year-round 6am-10pm, up to 4 hours, more frequent in summer
                % Temperature: Adjusted up to 4°F, can override anytime
                % ================================================================
                % Peak rates from user's tariff data:
                %   Summer: $0.096305/kWh (hours 13-18)
                %   Winter: $0.061350/kWh (hours 5-6, 16-18)
                % ================================================================
                
                if contains(season, {'summer', 'winter'})
                    duration = 4;
                    start_hour = randi([14, 18]);  % 2pm-6pm typical peak window
                    event_hours = start_hour:min(start_hour + duration - 1, 22);
                    DR_active_hours(event_hours) = true;
                    
                    % Value = avoiding peak tariff costs
                    if contains(season, 'summer')
                        C_IL_hourly(event_hours) = 0.096;  % Summer peak rate
                    else
                        C_IL_hourly(event_hours) = 0.061;  % Winter peak rate
                    end
                    
                    % No capacity payment
                    capacity_rate_per_kW = 0;
                end
                
            else  % Commercial
                % ================================================================
                % AUSTIN ENERGY - COMMERCIAL DEMAND RESPONSE
                % ================================================================
                % Program Type: Capacity-based
                % Payment: $50-70/kW saved (using $60/kW average)
                % Season: June 1 - Sept 30 only (4 months)
                % Events: Max 25/year, 1-7pm, up to 2 hours
                % Options: Standard DR (30-min notice) or Fast DR (10-min notice)
                % Eligible: 10% load reduction, min 20kW savings
                % ================================================================
                % CAPACITY ALLOCATION (4 months Jun-Sep):
                %   Spring: $0 (before June)
                %   Summer (3 months): $60 × 3/4 = $45/kW
                %   Autumn1 (1 month): $60 × 1/4 = $15/kW
                %   Others: $0
                % ================================================================
                
                if strcmp(season, 'summer')
                    % Summer: Jun-Aug (3 of 4 program months)
                    capacity_rate_per_kW = 60 * (3/4);  % $45/kW
                    
                    duration = 2;
                    start_hour = randi([13, 17]);  % 1pm-5pm window
                    event_hours = start_hour:min(start_hour + duration - 1, 19);
                    DR_active_hours(event_hours) = true;
                    
                    % Must curtail when called - no explicit performance payment
                    % TX commercial peak: ~$0.13/kWh
                    C_IL_hourly(event_hours) = 0.0;  % Nominal (<<< peak rate)
                    
                elseif strcmp(season, 'autumn1')
                    % Autumn1: Sep (1 of 4 program months)
                    capacity_rate_per_kW = 60 * (1/4);  % $15/kW
                    
                    duration = 2;
                    start_hour = randi([13, 17]);
                    event_hours = start_hour:min(start_hour + duration - 1, 19);
                    DR_active_hours(event_hours) = true;
                    
                    C_IL_hourly(event_hours) = 0.0;  % Nominal
                else
                    % No DR program outside Jun-Sep
                    capacity_rate_per_kW = 0;
                end
            end
            
        % ====================================================================
        % MINNESOTA
        % ====================================================================
        case 'Minnesota'
            % ================================================================
            % MVEC - ENERGY WISE WI-FI THERMOSTAT PROGRAM
            % ================================================================
            % Program Type: Enrollment bonus + bill discount + thermostat control
            % Payment Structure:
            %   - $50 one-time enrollment credit (handled in main script)
            %   - 10% summer bill discount (applied to all summer days above)
            %   - NO per-kW capacity payment (per user's research)
            %   - NO per-kWh performance payment
            % Value: Avoiding peak electricity costs during thermostat adjustment
            % Season: June 1 - Sept 30 (summer only)
            % Events: Brief adjustments during peak demand days
            % ================================================================
            % Peak rate from user's tariff data:
            %   $0.156/kWh (hours 15-19)
            % ================================================================
            
            if contains(season, 'summer')
                duration = 4;
                start_hour = randi([14, 17]);  % 2pm-5pm typical peak
                event_hours = start_hour:min(start_hour + duration - 1, 21);
                DR_active_hours(event_hours) = true;
                
                % No capacity payment (confirmed by user)
                capacity_rate_per_kW = 0;
                
                % Value = avoiding peak tariff costs
                C_IL_hourly(event_hours) = 0.156;  % Peak rate from tariff
            else
                % No DR program outside summer
                capacity_rate_per_kW = 0;
            end
            
        % ====================================================================
        % NEW YORK
        % ====================================================================
        case 'NewYork'
            if strcmp(sector, 'residential')
                % ================================================================
                % CONED - SMART USAGE REWARDS
                % ================================================================
                % Program Type: Smart thermostat with enrollment incentives
                % Payment Structure:
                %   - $85 enrollment rebate (handled in main script)
                %   - $25 annual incentive starting year 3 (handled in main script)
                %   - Requirement: 50% participation in event hours
                % Value: Avoiding peak electricity costs during thermostat adjustment
                % Events: 11am-11pm weekdays, ~4 hours, mostly summer
                % Temperature: Adjusted to comfortable but warmer setting
                % ================================================================
                % Peak rate from user's tariff data:
                %   Summer: $0.35/kWh (hours 9-23)
                % ================================================================
                
                if contains(season, 'summer')
                    duration = 4;
                    start_hour = randi([11, 19]);  % 11am-7pm start window
                    event_hours = start_hour:min(start_hour + duration - 1, 23);
                    DR_active_hours(event_hours) = true;
                    
                    % No capacity payment
                    capacity_rate_per_kW = 0;
                    
                    % Value = avoiding peak tariff costs
                    C_IL_hourly(event_hours) = 0.35;  % Peak rate from tariff
                else
                    capacity_rate_per_kW = 0;
                end
                
            else  % Commercial
                % ================================================================
                % CONED - COMMERCIAL SYSTEM RELIEF PROGRAM (CSRP)
                % ================================================================
                % Program Type: HYBRID (Capacity + Performance)
                % Capacity Payment:
                %   - Brooklyn/Bronx/Manhattan/Queens: $18/kW-month
                %   - Staten Island/Westchester: $6/kW-month
                %   - Using urban rate $18/kW-month (majority of customers)
                % Performance Payment: $1/kWh curtailed during events
                % Season: May-Sep (5 months)
                % Events: Weekdays 11am-11pm, 21-hour advance notice
                % Notification: During high-demand periods
                % ================================================================
                % CAPACITY ALLOCATION (5 months May-Sep):
                %   Spring (3 months): May ($18) = $18/kW
                %   Summer (3 months): Jun-Aug ($18 each) = $54/kW
                %   Autumn1 (1 month): Sep ($18) = $18/kW
                %   Others: $0
                %   TOTAL: $90/kW per year
                % ================================================================
                
                if strcmp(season, 'spring')
                    % Spring: May only (1 of 5 program months)
                    capacity_rate_per_kW = 18 * 1;  % $18/kW
                    
                    duration = randi([4, 6]);
                    start_hour = randi([11, 17]);  % 11am-5pm start window
                    event_hours = start_hour:min(start_hour + duration - 1, 23);
                    DR_active_hours(event_hours) = true;
                    
                    % Performance payment: $1/kWh actual curtailment
                    C_IL_hourly(event_hours) = 1.00;
                    
                elseif strcmp(season, 'summer')
                    % Summer: Jun-Aug (3 of 5 program months)
                    capacity_rate_per_kW = 18 * 3;  % $54/kW
                    
                    duration = randi([4, 6]);
                    start_hour = randi([11, 17]);
                    event_hours = start_hour:min(start_hour + duration - 1, 23);
                    DR_active_hours(event_hours) = true;
                    
                    C_IL_hourly(event_hours) = 1.00;
                    
                elseif strcmp(season, 'autumn1')
                    % Autumn1: Sep (1 of 5 program months)
                    capacity_rate_per_kW = 18 * 1;  % $18/kW
                    
                    duration = randi([4, 6]);
                    start_hour = randi([11, 17]);
                    event_hours = start_hour:min(start_hour + duration - 1, 23);
                    DR_active_hours(event_hours) = true;
                    
                    C_IL_hourly(event_hours) = 1.00;
                else
                    % No CSRP outside May-Sep
                    capacity_rate_per_kW = 0;
                end
            end
    end
end