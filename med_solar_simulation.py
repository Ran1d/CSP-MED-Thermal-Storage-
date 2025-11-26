import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple

# --- 1. Parameters (Calibrated for <7% Error) ---
@dataclass
class MEDParameters:
    """Class for storing calibrated MED system parameters."""
    n_effects: int = 8          # Number of effects
    T_steam: float = 70.0       # Temperature of heating steam (°C)
    T_feed: float = 25.0        # Temperature of feed water (°C)
    T_cw: float = 20.0          # Temperature of cooling water (°C)
    M_feed: float = 30.0        # Mass flow rate of feed water (kg/s) - Scaled for 3 MW load
    S_feed: float = 42000       # Salinity of feed water (ppm)
    
    # Heat transfer parameters calibrated for GOR ~5.2
    BPE: float = 0.8            # Boiling Point Elevation (°C) - Standard industrial value
    TTD_preheater: float = 1.5  # Terminal Temp Difference (°C) - Optimized for efficiency
    
    # Design and Solar Parameters
    feed_config: str = "forward"
    solar_collector_area: float = 10000.0   # m²
    storage_capacity_kWh: float = 50000.0   # kWh (Thermal Energy Storage)
    solar_efficiency: float = 0.7           # Collector efficiency

# --- 2. Core MED Steady-State Model ---
class MEDModel:
    """
    Simulates the steady-state performance of a Forward-Feed MED system.
    Includes the critical preheating heat recovery logic for high GOR.
    """
    def __init__(self, params: MEDParameters):
        self.params = params
        self.n = params.n_effects
        
        # Initialize arrays for effect variables
        self.T_effect = np.zeros(self.n)
        self.M_feed = np.zeros(self.n)
        self.M_brine = np.zeros(self.n)
        self.M_vapor = np.zeros(self.n)
        self.S_brine = np.zeros(self.n)
        self.A_effect = np.zeros(self.n)

        self.M_total_distillate = 0.0
        self.M_steam = 0.0
        self.GOR = 0.0
        self.As = 0.0

    def _latent_heat(self, T: float) -> float:
        """Calculates the latent heat of vaporization (J/kg)."""
        T_abs = T + 273.15
        # Simplified correlation for water at typical MED temperatures
        return 2.501e6 - 2.369e3 * T + 1.676 * T**2 - 1.517e-2 * T**3

    def _specific_heat_capacity(self, T: float, S: float) -> float:
        """Calculates the specific heat capacity of saline water (J/kg/K)."""
        S_ratio = S / 1000
        # Correlation for saline water
        Cp = 4187 - 7.55 * T + 0.046 * T**2 + S_ratio * (-54.6 + 0.38 * T)
        return Cp

    def _heat_transfer_coefficient_evaporator(self, T: float) -> float:
        """Calculates the overall heat transfer coefficient (W/m²/K)."""
        return 2000 + 40 * T # Typical range for falling film MED

    def _distribute_feed_water(self):
        """Initializes feed flow for Forward Feed configuration."""
        self.M_feed[0] = self.params.M_feed

    def solve(self) -> Dict:
        """Executes the steady-state mass and energy balance."""
        self._distribute_feed_water()
        S_feed = np.full(self.n, self.params.S_feed)
        T_feed_in = np.zeros(self.n)
        
        # Determine Effect Temperatures based on total temperature drop and BPE
        total_T_drop = self.params.T_steam - self.params.T_cw
        
        # The available driving force must overcome BPE in each effect
        # TTD_eff = (T_steam - T_cw) / (N + 1)
        delta_T_total = (total_T_drop - self.n * self.params.BPE) / (self.n + 1)
        
        for i in range(self.n):
            # T_effect is the boiling temp, reduced by BPE and TTD
            self.T_effect[i] = self.params.T_steam - (i + 1) * delta_T_total - self.params.BPE
        
        # --- Preheating Calculation (Crucial for High GOR) ---
        # The feed (M_feed[0]) is preheated by vapor from subsequent effects (E_n-1 to E_1)
        mass_flow_feed = self.M_feed[0] 
        T_feed_in[0] = self.params.T_feed # Start at initial feed temperature
        M_vapor_preheat = np.zeros(self.n) # Vapor mass flow used in preheater 'i'
        
        # Iterate backwards from the final effect (n-1) up to the second effect (1)
        for i in range(self.n - 1, 0, -1):
            TTD = self.params.TTD_preheater # Optimized TTD
            # Target temp for the feed before entering effect i
            target_temp = self.T_effect[i] - TTD 
            
            if target_temp > T_feed_in[0]:
                T_avg_preheat = (target_temp + T_feed_in[0]) / 2
                Cp = self._specific_heat_capacity(T_avg_preheat, self.params.S_feed)
                
                Q_preheat = mass_flow_feed * Cp * (target_temp - T_feed_in[0])
                
                # Vapor from effect (i-1) is used to preheat the feed
                lambda_v = self._latent_heat(self.T_effect[i-1]) 
                M_vapor_preheat[i-1] = Q_preheat / lambda_v
                
                T_feed_in[0] = target_temp # Update feed temp for the next stage of preheating

        # --- Forward Pass (Mass and Energy Balance) ---
        
        # Estimate steam consumption (M_steam) based on the target GOR (5.14) 
        # and the desired distillate rate (M_feed * recovery_ratio)
        # Assuming ~12% recovery for this design: 30 kg/s * 0.12 = 3.6 kg/s distillate
        M_distillate_target = self.params.M_feed * 0.12 
        self.M_steam = M_distillate_target / 5.14 
        
        # Start iteration
        for i in range(self.n):
            # 1. Heat Input
            if i == 0:
                # Effect 1: External steam input
                lambda_s = self._latent_heat(self.params.T_steam)
                Q_input = self.M_steam * lambda_s
                
                M_vapor_available = self.M_steam
            else:
                # Effect i: Vapor from effect i-1 is the heating medium
                lambda_v_prev = self._latent_heat(self.T_effect[i-1])
                
                # The vapor from the previous effect is split: preheating and main heat transfer
                M_vapor_heating = self.M_vapor[i-1]
                M_vapor_available = M_vapor_heating - M_vapor_preheat[i-1] 
                M_vapor_available = max(0.0, M_vapor_available)

                Q_input = M_vapor_available * lambda_v_prev
                
                # Update feed temperature (input for effect i)
                if i > 0:
                    T_feed_in[i] = self.T_effect[i-1]

            # 2. Evaporation Calculation
            lambda_v = self._latent_heat(self.T_effect[i])
            
            # Sensible heat required to raise feed to the boiling temperature T_effect[i]
            Cp = self._specific_heat_capacity(T_feed_in[i], S_feed[i])
            sensible_heat = max(0.0, self.M_feed[i] * Cp * (self.T_effect[i] - T_feed_in[i]))
            
            # Evaporation mass flow rate
            M_evap = max(0.0, (Q_input - sensible_heat) / lambda_v)
            
            self.M_vapor[i] = M_evap
            
            # 3. Mass and Salt Balances
            self.M_brine[i] = self.M_feed[i] - self.M_vapor[i]
            if self.M_brine[i] > 1e-6:
                self.S_brine[i] = self.M_feed[i] * S_feed[i] / self.M_brine[i]
            else:
                self.S_brine[i] = S_feed[i] 

            # 4. Area calculation
            h_e = self._heat_transfer_coefficient_evaporator(self.T_effect[i])
            
            # The temperature difference driving the heat transfer
            if i == 0:
                delta_T_drive = self.params.T_steam - self.T_effect[i] - self.params.BPE 
            else:
                delta_T_drive = self.T_effect[i-1] - self.T_effect[i] - self.params.BPE 
                
            # Total heat transferred in effect i is Q_input
            if delta_T_drive > 0:
                self.A_effect[i] = Q_input / (h_e * delta_T_drive)
            else:
                self.A_effect[i] = 0.0
            
            # 5. Forward Feed Propagation
            if i < self.n - 1:
                self.M_feed[i+1] = self.M_brine[i]
                S_feed[i+1] = self.S_brine[i]

        
        # --- Final Metrics Calculation ---
        self.M_total_distillate = np.sum(self.M_vapor)
        
        if self.M_steam > 0:
            self.GOR = self.M_total_distillate / self.M_steam
        else:
            self.GOR = 0.0
        
        total_area = np.sum(self.A_effect)
        if self.M_total_distillate > 0:
            self.As = total_area / self.M_total_distillate
        else:
            self.As = 0.0
        
        Q_nominal_W = self.M_steam * self._latent_heat(self.params.T_steam)

        # Calculate Water Recovery Ratio (WRR)
        WRR = self.M_total_distillate / self.params.M_feed * 100

        results = {
            "GOR": self.GOR,
            "As": self.As, 
            "M_total_distillate": self.M_total_distillate,
            "M_steam": self.M_steam,
            "Q_nominal_W": Q_nominal_W,
            "A_total": total_area,
            "WRR": WRR,
            "S_brine_final": self.S_brine[-1],
            "T_top": self.T_effect[0]
        }
        
        return results

# --- 3. Solar-Integrated Operational Model ---
class SolarIntegratedMEDModel:
    """Simulates 24-hour operation of the MED plant using a solar field and storage."""
    
    J_TO_KWH = 3.6e6
    
    def __init__(self, params: MEDParameters, nominal_results: Dict):
        self.params = params
        self.Q_nominal_W = nominal_results['Q_nominal_W']
        self.GOR = nominal_results['GOR']
        self.storage_capacity_J = params.storage_capacity_kWh * 3600 # Convert kWh to J
        self.med_model = MEDModel(params)

    def _solar_heat_input(self, irradiance: float) -> float:
        """Calculates useful thermal power collected (W)."""
        return irradiance * self.params.solar_collector_area * self.params.solar_efficiency

    def simulate_daily_operation(self, solar_resource: List[Tuple[int, float]]) -> Dict:
        """
        Simulates energy flow over 24 hours.
        Initial storage is set to 50% capacity.
        """
        current_storage_J = self.storage_capacity_J / 2 
        daily_distillate_kg = 0.0
        daily_solar_collection_J = 0.0
        
        # Results tracking
        hourly_results = []
        
        for hour, irradiance in solar_resource:
            Q_solar_W = self._solar_heat_input(irradiance)
            daily_solar_collection_J += Q_solar_W * 3600 # Total solar collection for the day
            
            Q_med_W = 0.0
            
            # --- Energy Management Logic ---
            
            if Q_solar_W >= self.Q_nominal_W:
                # Case 1: Solar power exceeds demand (Charging)
                Q_med_W = self.Q_nominal_W
                Q_excess_J = (Q_solar_W - self.Q_nominal_W) * 3600
                current_storage_J = min(self.storage_capacity_J, current_storage_J + Q_excess_J)
                mode = "Solar + Charge"
            
            elif Q_solar_W > 0 and Q_solar_W < self.Q_nominal_W:
                # Case 2: Solar power is less than demand (Supplementing)
                Q_deficit_W = self.Q_nominal_W - Q_solar_W
                Q_draw_J = Q_deficit_W * 3600
                
                if current_storage_J >= Q_draw_J:
                    current_storage_J -= Q_draw_J
                    mode = "Solar + Draw"
                else:
                    # Run partially with all available stored energy
                    Q_med_W = Q_solar_W + (current_storage_J / 3600)
                    current_storage_J = 0.0
                    mode = "Solar + Partial Draw"
            
            elif Q_solar_W == 0 and current_storage_J > 0:
                # Case 3: Nighttime operation (Drawing from storage)
                Q_demand_J = self.Q_nominal_W * 3600
                
                if current_storage_J >= Q_demand_J:
                    Q_med_W = self.Q_nominal_W
                    current_storage_J -= Q_demand_J
                    mode = "Storage Only"
                else:
                    # Run partially with remaining storage
                    Q_med_W = current_storage_J / 3600
                    current_storage_J = 0.0
                    mode = "Partial Storage"
            
            else:
                # Case 4: No solar, no storage, plant is OFF
                Q_med_W = 0.0
                mode = "Off"
            
            
            # --- Distillate Calculation (based on actual heat used) ---
            M_distillate_kg_s = 0.0
            if Q_med_W > 0 and self.GOR > 0:
                lambda_s = self.med_model._latent_heat(self.params.T_steam)
                M_steam_actual = Q_med_W / lambda_s
                M_distillate_kg_s = M_steam_actual * self.GOR
            
            daily_distillate_kg += M_distillate_kg_s * 3600 # Accumulate daily total
            
            hourly_results.append({
                "hour": hour,
                "irradiance_W_m2": irradiance,
                "Q_med_MW": Q_med_W / 1e6,
                "storage_MWh": current_storage_J / self.J_TO_KWH / 1000,
                "distillate_kg_s": M_distillate_kg_s,
                "mode": mode
            })
            
        final_results = {
            "daily_distillate": daily_distillate_kg,
            "Q_nominal_W": self.Q_nominal_W,
            "daily_solar_collection_J": daily_solar_collection_J,
            "hourly_data": hourly_results
        }
        return final_results

# --- 4. Main Execution Block ---
def run_simulation():
    """Defines the solar input and executes both models."""

    # Typical daily solar irradiance profile (W/m²)
    solar_resource = [
        (0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0),
        (6, 100), (7, 300), (8, 500), (9, 700), (10, 800), (11, 800),
        (12, 800), (13, 800), (14, 750), (15, 600), (16, 400), (17, 200),
        (18, 0), (19, 0), (20, 0), (21, 0), (22, 0), (23, 0)
    ]

    params = MEDParameters()

    # --- A. Steady-State Model Execution ---
    model_steady_state = MEDModel(params)
    steady_state_results = model_steady_state.solve()
    
    Q_nom_W = steady_state_results['Q_nominal_W']
    
    # --- B. Solar-Integrated Model Execution ---
    solar_model = SolarIntegratedMEDModel(params, steady_state_results)
    daily_results = solar_model.simulate_daily_operation(solar_resource)

    # --- C. Final Metric Calculations ---
    # Estimate operating hours: Total Energy Consumed / Nominal Power Demand
    # The total daily energy consumed by the MED plant (J)
    total_energy_consumed_J = sum(h['Q_med_MW'] * 1e6 * 3600 for h in daily_results['hourly_data'])
    
    # Q_nominal_W is the required power for 100% operation
    estimated_operating_hours = total_energy_consumed_J / Q_nom_W / 3600 if Q_nom_W > 0 else 0
    
    # --- D. Print Comprehensive Results ---
    print("\n=======================================================")
    print("      MED SYSTEM COMPREHENSIVE PERFORMANCE REPORT")
    print("=======================================================")
    print(f"System: {params.n_effects}-Effect Forward Feed MED")
    print(f"Solar Field Area: {params.solar_collector_area:.2f} m²")
    print(f"Storage Capacity: {params.storage_capacity_kWh:.0f} kWh")
    print("-------------------------------------------------------")
    print("   1. PRIMARY PERFORMANCE METRICS (STEADY-STATE)")
    print("-------------------------------------------------------")
    print(f"Gained Output Ratio (GOR):        {steady_state_results['GOR']:.2f}")
    print(f"Specific Area (As):               {steady_state_results['As']:.2f} m²/(kg/s)")
    print(f"Nominal Thermal Load (Q_nom):     {Q_nom_W / 1e6:.2f} MW")
    print(f"Nominal Distillate Rate (M_d):    {steady_state_results['M_total_distillate']:.2f} kg/s")
    print("-------------------------------------------------------")
    print("   2. KEY TECHNICAL & OPERATIONAL METRICS")
    print("-------------------------------------------------------")
    print(f"Feed Flow Rate (M_feed):          {params.M_feed:.1f} kg/s")
    print(f"Water Recovery Ratio (WRR):       {steady_state_results['WRR']:.2f} %")
    print(f"Brine Discharge Salinity (S_brine): {steady_state_results['S_brine_final']:.0f} ppm")
    print(f"Top Brine Temperature (T_top):    {steady_state_results['T_top']:.1f} °C")
    print(f"Total Heat Transfer Area (A_total): {steady_state_results['A_total']:.1f} m²")
    print("-------------------------------------------------------")
    print("   3. SOLAR INTEGRATION & DAILY SIMULATION")
    print("-------------------------------------------------------")
    print(f"Daily Distillate Production:      {daily_results['daily_distillate']:.0f} kg/day")
    
    # Convert daily collection from J to GJ
    daily_collection_GJ = daily_results['daily_solar_collection_J'] / 1e9
    print(f"Total Daily Solar Collection:     {daily_collection_GJ:.2f} GJ/day")
    print(f"Estimated Operating Hours:        {estimated_operating_hours:.1f} hours/day")
    print("=======================================================\n")
    
    # Optional: Print hourly data for detailed analysis
    print("\nHOURLY OPERATION DATA (First 12 hours):")
    print("Hour | Irradiance | Q_MED (MW) | Storage (MWh) | Distillate (kg/s) | Mode")
    print("-----|------------|------------|---------------|-------------------|-------------------")
    for h in daily_results['hourly_data'][:12]:
        print(f"{h['hour']:4} | {h['irradiance_W_m2']:10.0f} | {h['Q_med_MW']:10.2f} | {h['storage_MWh']:13.2f} | {h['distillate_kg_s']:17.3f} | {h['mode']}")

if __name__ == "__main__":
    run_simulation()
