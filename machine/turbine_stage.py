import math
from dataclasses import dataclass

from configs.constants import *
from configs.modes import *
from core.geometry_models import SectionDiameters
from core.gas_func import velocity_critical, pi_lambda, q_lambda

@dataclass
class TurbineStageParameters:
    mode: str

    G_in: float
    G_out: float
    h_2_rotor: float
    u_mid: float
    inlet: SectionDiameters
    D_hub_out: float

    p_in: float
    pi: float
    T_in: float
    T_out:float
    L_st:float

    reaction: float
    y: float
    S_stator: float
    S_rotor: float

    r_stator: float
    r_rotor: float
    r_stator_colling: float | None = None
    r_rotor_colling: float | None = None

    phi_cooling: float | None = None
    psi_cooling: float | None = None

@dataclass
class BladeRowGridParams:
    t_opt_rel: float
    delta_t_opt_rel: float
    gamma_deg: float
    blade_chord: float
    r_out: float
    S_out_rel: float
    K_crit: float
    solidity: float
    blade_count: int

class TurbineStage:
    def __init__(self, params:TurbineStageParameters):
        self.params = params

    def calculate(self):
        h_0_stator = (self.params.inlet.tip - self.params.inlet.hub) / 2 
        h_1_stator = h_0_stator + (self.params.h_2_rotor - h_0_stator) / 2

        L_ad = 0.5 * (self.params.u_mid / self.params.y) ** 2

        c_ad = 1.415 * math.sqrt(L_ad)

        a_1crit = velocity_critical(k=K_GAS, R=R_GAS, T=self.params.T_in)
        lambda_ad = c_ad / a_1crit

        pi_lambda_ad = pi_lambda(lam=lambda_ad, k=K_GAS)

        p_2_stat = pi_lambda_ad * self.params.p_in

        y = self.params.u_mid / c_ad

        L_ad_stator = (1 - self.params.reaction) * L_ad
        L_ad_rotor = self.params.reaction * L_ad

        c_1t = math.sqrt(2 * L_ad_stator)
        c_1 = self.params.phi_cooling * c_1t

        T_1 = self.params.T_in - self.params.phi_cooling ** 2 * L_ad_stator / (K_GAS / (K_GAS - 1) * R_GAS)

        lambda_1t = c_1t / a_1crit
        pi_lambda_1t = pi_lambda(lam=lambda_1t, k=K_GAS)
        p_1_stat = pi_lambda_1t * self.params.p_in
        density_1 = p_1_stat / (R_GAS * T_1)

        alpha_1_rad = math.asin(self.params.G_in / (math.pi * self.params.inlet.mid * h_1_stator * c_1 * density_1))
        alpha_1_deg = math.degrees(alpha_1_rad)

        K_reaction = ((self.params.inlet.mid / (self.params.inlet.mid - h_1_stator)) * math.cos(alpha_1_rad)) ** 2 + math.sin(alpha_1_rad) ** 2
        reaction_hub = 1 - (1 - self.params.reaction) * K_reaction

        w_1 = math.sqrt(c_1**2 + self.params.u_mid**2 - 2 * c_1 * self.params.u_mid * math.cos(alpha_1_rad))
        beta_1_rad = math.asin(c_1 * math.sin(alpha_1_rad) / w_1)
        beta_1_deg = math.degrees(beta_1_rad)

        w_2 = self.params.psi_cooling * math.sqrt(w_1**2 + 2 * L_ad_rotor)

        T_w1 = self.params.T_in - (c_1**2 - w_1**2) / (2 * (K_GAS / (K_GAS - 1) * R_GAS))

        a_w1crit = velocity_critical(k=K_GAS, R=R_GAS, T=T_w1)
        lambda_w2 = w_2 / a_w1crit
        pi_lambda_w2 = pi_lambda(lam=lambda_w2, k=K_GAS)
        q_lambda_w2 = q_lambda(lam=lambda_w2, k=K_GAS)
        p_w2 = p_2_stat / pi_lambda_w2

        beta_2_rad = math.asin(self.params.G_out * math.sqrt(T_w1) / 
                               (math.pi * self.params.inlet.mid * self.params.h_2_rotor * q_lambda_w2 * S_GAS * p_w2))
        beta_2_deg = math.degrees(beta_2_rad)

        c_2 = math.sqrt(w_2**2 + self.params.u_mid**2 - 2 * w_2 * self.params.u_mid * math.cos(beta_2_rad))

        alpha_2_rad = math.asin(w_2 * math.sin(beta_2_rad) / c_2)
        alpha_2_deg = math.degrees(alpha_2_rad)
        T_2 = self.params.T_in -  self.params.L_st / (K_GAS / (K_GAS - 1) * R_GAS)
        a_2crit = velocity_critical(k=K_GAS, R=R_GAS, T=T_2)
        lambda_c2 = c_2 / a_2crit
        pi_lambda_c2 = pi_lambda(lam=lambda_c2, k=K_GAS)

        p_2 = p_2_stat / pi_lambda_c2

        delta_beta = 180 - (beta_1_deg + beta_2_deg)
        delta_alpha = 180 - (90 + alpha_1_deg)

        K_stator = 1 / math.sin(alpha_1_rad)
        K_rotor = math.sin(beta_1_rad) / math.sin(beta_2_rad)

        stator_grid = self.blade_row_grid_params(
            K=K_stator,
            delta_angle_deg=delta_alpha,
            lambda_value=lambda_1t,
            gamma_angle_deg=90 - alpha_1_deg,
            S=self.params.S_stator,
            r=self.params.r_stator,
            r_cooling=self.params.r_stator_colling
        )

        rotor_grid = self.blade_row_grid_params(
            K=K_rotor,
            delta_angle_deg=delta_beta,
            lambda_value=lambda_w2 / self.params.psi_cooling,
            gamma_angle_deg=beta_1_deg - beta_2_deg,
            S=self.params.S_rotor,
            r=self.params.r_rotor,
            r_cooling=self.params.r_rotor_colling
        )

        solidity_hub_rotor = math.pi * self.params.inlet.hub / rotor_grid.blade_count


        print(f'h_1 {h_1_stator}')
        print(f'Work ad {L_ad}')
        print(f'lambda ad {lambda_ad}')
        print(f'pressure out rotor {p_2_stat}')
        print(f'c 1 {c_1}')
        print(f'temp out stator {T_1}')
        print(f'alpha 1 {alpha_1_deg}')
        print(f'K = {K_reaction}')
        print(f'w_1 = {w_1}')
        print(f'beta 1 {beta_1_deg}')
        print(f'w2 = {w_2}')
        print(f'Tw1 = {T_w1}')
        print(f'beta 2 {beta_2_deg}')
        print(f'c 2 = {c_2}')
        print(f'alpha 2 {alpha_2_deg}')
        print(f'T out {T_2, self.params.T_out}')
        print(f'lambda c2 {lambda_c2}')
        print(f'p out {p_2}')
        print(f'delta beta {delta_beta}')
        print(f'delta alpha {delta_alpha}')
        print(f'K stator, rotor {K_stator, K_rotor}')
        print(f'решетка статора {stator_grid}')
        print(f'решетка ротора {rotor_grid}')



    def blade_row_grid_params(self, K: float,delta_angle_deg: float,lambda_value: float,gamma_angle_deg: float,
        S: float,r: float, r_cooling:float) -> BladeRowGridParams:
        
        def optimal_grid_step(K: float, delta_deg: float) -> float:
            delta_rad = math.radians(delta_deg)
            if K >= 1.5:
                return (
                    0.327 / (K ** 0.327 * math.cbrt(delta_rad))
                    - 0.994 / K ** 0.385
                    + 1.314)
            if K > 1:
                return (
                    (1.727 / K - 0.869) / math.cbrt(delta_rad)
                    - 1.71 / K
                    + 1.604)
            raise ValueError(f"K must be > 1, got K = {K}")

        def optimal_grid_step_delta(lam: float) -> float:
            return -0.625 * lam ** 2 + 0.48 * lam + 0.016

        def gamma_calc(angle_deg: float) -> float:
            return 68.7 + 9.33e-4 * angle_deg - 6.052e-3 * angle_deg ** 2

        def trailing_edge_thickness_coeff(S_rel: float, t_rel: float) -> float:
            return 1 - 15 * S_rel ** 2 + (3.75 * t_rel - 0.6) * S_rel
        
        t_opt_rel = optimal_grid_step(K=K, delta_deg=delta_angle_deg)
        delta_t_opt_rel = optimal_grid_step_delta(lam=lambda_value)

        gamma_deg = gamma_calc(angle_deg=gamma_angle_deg)

        blade_chord = S / math.sin(math.radians(gamma_deg))
        r_out = r * blade_chord

        S_out_rel = 2 * r_out / blade_chord

        K_crit = trailing_edge_thickness_coeff(S_rel=S_out_rel,t_rel=t_opt_rel)

        solidity = (1 + delta_t_opt_rel) * K_crit * t_opt_rel * blade_chord

        solidity_cooling = solidity * r_cooling

        z = round(math.pi * self.params.inlet.mid / solidity_cooling)

        return BladeRowGridParams(
            t_opt_rel=t_opt_rel,
            delta_t_opt_rel=delta_t_opt_rel,
            gamma_deg=gamma_deg,
            blade_chord=blade_chord,
            r_out=r_out,
            S_out_rel=S_out_rel,
            K_crit=K_crit,
            solidity=solidity_cooling,
            blade_count=z
            )

