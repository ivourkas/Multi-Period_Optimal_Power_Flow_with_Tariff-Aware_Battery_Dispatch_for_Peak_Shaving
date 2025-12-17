function [res_event_days, com_event_days, overlap_days, total_weekdays] = ...
    get_DR_event_frequency(location, season, is_residential)
    % ========================================================================
    % DR EVENT FREQUENCY CALCULATOR - INDEPENDENT SECTOR EVENTS
    % ========================================================================
    % Returns event frequencies for EACH SECTOR INDEPENDENTLY
    % Recognizes that residential and commercial DR events occur on different days
    %
    % INPUTS:
    %   location: 'California', 'Texas', 'Minnesota', 'NewYork'
    %   season: 'summer', 'winter', 'spring', 'autumn1', 'autumn2'
    %   is_residential: Boolean array (1×10) - used to check if sector exists
    %
    % OUTPUTS:
    %   res_event_days: Number of weekdays with residential DR events
    %   com_event_days: Number of weekdays with commercial DR events
    %   overlap_days: Days where BOTH sectors have events (usually 0-2)
    %   total_weekdays: Total weekdays in season
    %
    % NOTE: Actual scenario days calculated as:
    %   - Res only: res_event_days - overlap_days
    %   - Com only: com_event_days - overlap_days
    %   - Both: overlap_days
    %   - None: total_weekdays - res_event_days - com_event_days + overlap_days
    % ========================================================================
    
    % Calculate weekdays in season
    if contains(season, {'summer', 'winter', 'spring'})
        weeks = 12;
    elseif strcmp(season, 'autumn1')
        weeks = 4;
    else  % autumn2
        weeks = 8;
    end
    total_weekdays = 5 * weeks;
    
    % Check if sectors exist in community (excluding slack bus)
    has_residential = sum(is_residential(2:end)) > 0;
    has_commercial = sum(~is_residential(2:end)) > 0;
    
    % ========================================================================
    % GET EVENT FREQUENCIES FOR EACH SECTOR
    % ========================================================================
    if has_residential
        res_event_days = get_sector_frequency(location, season, true, weeks);
    else
        res_event_days = 0;
    end
    
    if has_commercial
        com_event_days = get_sector_frequency(location, season, false, weeks);
    else
        com_event_days = 0;
    end
    
    % ========================================================================
    % ESTIMATE OVERLAP (BOTH EVENTS SAME DAY)
    % ========================================================================
    % Most DR programs are called independently, but occasional overlap
    % Conservative estimate: ~10% overlap when both programs active
    % Exception: If one program has no events, overlap = 0
    
    if res_event_days > 0 && com_event_days > 0
        % Both programs have events - estimate small overlap
        % Use minimum frequency as upper bound for overlap
        max_possible_overlap = min(res_event_days, com_event_days);
        
        % Estimate 10% overlap (conservative)
        overlap_days = round(0.10 * max_possible_overlap);
        
        % Ensure overlap doesn't exceed either event frequency
        overlap_days = min(overlap_days, min(res_event_days, com_event_days));
    else
        % At least one program has no events
        overlap_days = 0;
    end
    
    % Ensure frequencies don't exceed total weekdays
    res_event_days = min(res_event_days, total_weekdays);
    com_event_days = min(com_event_days, total_weekdays);
    
    % Final sanity check: total events can't exceed weekdays
    if (res_event_days + com_event_days - overlap_days) > total_weekdays
        % Reduce overlap to make it fit
        overlap_days = max(0, res_event_days + com_event_days - total_weekdays);
    end
end

% ============================================================================
% HELPER FUNCTION: GET FREQUENCY FOR SINGLE SECTOR
% ============================================================================
function event_days = get_sector_frequency(location, season, is_residential_sector, weeks)
    % Returns event frequency for a single sector
    
    switch location
        case 'California'
            if is_residential_sector
                % Power Saver Rewards: 12 events/year, May 1 - Oct 31
                % Program active months: May, Jun, Jul, Aug, Sep, Oct (6 months = 24 weeks)
                % But our seasons don't align perfectly with calendar months
                
                if strcmp(season, 'spring')
                    % Spring = Mar-May (12 weeks)
                    % Only May (last 1/3) is in program ≈ 4 weeks
                    event_days = round(12 * 4 / 24);  % = 2 events
                    
                elseif strcmp(season, 'summer')
                    % Summer = Jun-Aug (12 weeks), fully in program
                    event_days = round(12 * 12 / 24);  % = 6 events
                    
                elseif strcmp(season, 'autumn1')
                    % Autumn1 = Sep (4 weeks), fully in program
                    event_days = round(12 * 4 / 24);  % = 2 events
                    
                elseif strcmp(season, 'autumn2')
                    % Autumn2 = Oct-Nov (8 weeks)
                    % Only Oct (first 1/2) is in program ≈ 4 weeks
                    event_days = round(12 * 4 / 24);  % = 2 events
                    
                else  % winter
                    event_days = 0;
                end
                
            else  % Commercial BIP
                % BIP: Year-round, 30 events/year (180 hours ÷ 6 hours/event)
                % Events concentrated during high grid stress periods
                
                if strcmp(season, 'summer')
                    % Summer: Peak demand, heat waves, flex alerts
                    event_days = 12;  % 40% of annual events
                    
                elseif strcmp(season, 'spring')
                    % Spring: Shoulder season, moderate demand
                    event_days = 6;  % 20% of annual events
                    
                elseif strcmp(season, 'autumn1')
                    % September: Still hot in California, high AC use
                    event_days = 4;  % 13% of annual events
                    
                elseif strcmp(season, 'autumn2')
                    % October-November: Cooling down
                    event_days = 4;  % 13% of annual events
                    
                else  % winter
                    % Winter: Lowest stress period
                    event_days = 4;  % 13% of annual events
                end
                
                % Total: 12 + 6 + 4 + 4 + 4 = 30 events ✓
            end
        case 'Texas'
            if is_residential_sector
                % Direct Energy: Year-round
                if strcmp(season, 'summer')
                    event_days = 15;
                elseif strcmp(season, 'winter')
                    event_days = 5;
                else
                    event_days = 0;
                end
            else  % Commercial
                % Austin Energy: Jun-Sep only
                if strcmp(season, 'summer')
                    event_days = round(25 * 0.75);  % ~19 events
                elseif strcmp(season, 'autumn1')
                    event_days = round(25 * 0.25);  % ~6 events
                else
                    event_days = 0;
                end
            end
            
        case 'Minnesota'
            % MVEC: Same for both sectors
            if strcmp(season, 'summer')
                event_days = 5;
            else
                event_days = 0;
            end
            
        case 'NewYork'
            if is_residential_sector
                % Smart Usage Rewards: Summer only
                if strcmp(season, 'summer')
                    event_days = 10;
                else
                    event_days = 0;
                end
            else  % Commercial
                % CSRP: May-Sep
                if strcmp(season, 'spring')
                    event_days = 3;
                elseif strcmp(season, 'summer')
                    event_days = 10;
                elseif strcmp(season, 'autumn1')
                    event_days = 3;
                else
                    event_days = 0;
                end
            end
            
        otherwise
            event_days = 0;
    end
end