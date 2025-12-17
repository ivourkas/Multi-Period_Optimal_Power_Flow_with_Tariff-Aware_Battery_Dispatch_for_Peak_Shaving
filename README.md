# Battery Storage Optimization with Demand Response

MATLAB framework for multi-year techno-economic analysis of distributed battery energy storage systems with demand response participation.

## Overview

This framework simulates 10-year operation of battery storage systems across a 10-bus distribution network, evaluating economic performance under various configurations of PV penetration, battery sizing, demand response participation, and customer mix (residential/commercial). Two degradation models are supported: simple throughput-based and rainflow cycle counting.

## Features

- **10-year simulation** with seasonal resolution (spring, summer, autumn, winter)
- **Battery degradation tracking** with automatic replacement scheduling
- **Rainflow cycle counting** for accurate Li-ion/Na-ion lifetime estimation
- **Demand response** event generation and revenue calculation
- **Location-specific** temperature and tariff profiles (CA, TX, NY, MN)
- **Constrained OPF** with battery dispatch optimization
- **Parallel execution** for parameter sweeps

## Repository Structure

```
├── simple_full_v2.m              # Main simulation (simple degradation)
├── rainflow_full_v2.m            # Main simulation (rainflow degradation)
├── COPF_fe.m                     # Constrained optimal power flow solver
├── BattPkShvRsCm_v2.m            # Battery peak shaving model
│
├── define_tariffs.m              # TOU tariff structures
├── get_tariff.m                  # Tariff lookup
├── get_battery_characteristics.m # Li-ion / Na-ion parameters
├── get_seasonal_temperature.m    # Location temperature profiles
│
├── generate_DR_events_all_locations.m  # DR event generator
└── get_DR_event_frequency.m            # DR event parameters
```

## Parameter Space

| Parameter | Values |
|-----------|--------|
| PV Scale | 0.2, 0.5, 1.2 |
| Battery Type | Li-ion, Na-ion |
| Battery Scale | 0.0, 0.5, 1.0, 1.5 |
| Residential Mix | 3/9, 6/9, 9/9 |
| IL Scale | 0.05, 0.1, 0.3 |
| Peak Penalty | 0.01, 0.1, 1, 10 |
| Location | California, Texas, New York, Minnesota |

## Quick Start

```matlab
% Run simple degradation model (Texas)
simple_full_v2

% Run rainflow degradation model (California)
rainflow_full_v2
```

## Outputs

- `simple_full_TX.csv` — Consolidated results for all parameter combinations
- `pk_per_bus_simple_full_TX.csv` — Per-bus worst-case peak values by season

**Metrics tracked:**
- Discounted objective function value (NPV)
- Total load served
- Battery replacements per bus
- DR events and revenue
- Seasonal peak demands

## Degradation Models

| Model | File | Method |
|-------|------|--------|
| Simple | `simple_full_v2.m` | Throughput-based cycle counting |
| Rainflow | `rainflow_full_v2.m` | ASTM E1049 rainflow algorithm |

## Requirements

- MATLAB R2020a or later
- Parallel Computing Toolbox (optional, for speedup)
- Optimization Toolbox

## License

MIT License
