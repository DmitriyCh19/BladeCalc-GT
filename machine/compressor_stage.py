import math
from dataclasses import dataclass

from configs.constants import *
from configs.modes import *
from core.gas_func import velocity_critical, q_lambda, local_speed_velocity
from core.geometry import area_calc, relative_diameter_hub, calculate_section_diameters
from core.geometry_models import VelocityTriangle, RotorGeometry


# D_tip = const
z = 9

alpha11_deg = 60.546

L_stage_rel = [0.24227, 0.252, 0.271, 0.281, 0.271, 0.261, 0.251, 0.241]
D_tip = 0.5923
eff_stage = [0.87, 0.89, 0.89, 0.9, 0.91, 0.9, 0.885, 0.885]
eff_rel = [1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.935, 0.93]
reaction_stage = [0.5, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57]
c1a_stage = [193.4903, 191.9, 190.3, 188.7, 187.2, 187.2, 185.6, 184.01, 182.4]
h_rot_rel = [3.5, 3.5, 3.3, 3.1, 2.9, 2.7, 2.4, 2.2, 2]
G = 72.4359
d_hub_rel = 0.74
lambda_in = 0.58
u_tip_1 = 393.4265
n = 12692
T_in = 437.9788
p_in = 355623
pi = 6.58
D_1mid_1 = 0.521
D_1hub_1 = 0.4383

@dataclass
class CompressorStageParameters:
    mode: str
    D_const: float

    alpha1_deg: float
    L_rel: float
    D_tip: float
    efficiency: float
    eff_tip_rel: float
    reaction: float
    c_1a: float
    c_3a: float
    h_rot_rel: float
    G: float
    d_hub_rel: float
    lambda_in: float
    u_tip: float
    T_in: float
    p_in: float
    D_1mid: float
    D_1hub: float
    K_g: float
    eff_rotor: float

@dataclass
class StageThermodynamics:
    L_stage: float

    T_in: float
    T_out: float

    p_in: float
    p_2:float
    p_out: float

    pi: float

@dataclass
class InletSolution:
    area: float

    u_mid: float

    lambda_1: float

    L_ku_rel: float

class CompressorStage:
    def __init__(self, params:CompressorStageParameters):
        self.params = params

    def calculate(self):
        thermo = self.calculate_thermodynamics()
        
        L_ku = thermo.L_stage* self.params.eff_tip_rel

        u_mid_in = self.params.u_tip * math.sqrt((1 + self.params.d_hub_rel ** 2) / 2)
        a_1crit = velocity_critical(k=K_AIR, R=R_AIR, T=thermo.T_in)

        # пропущенны пункты 9-10

        lambda_1a = self.params.c_1a / a_1crit
        q_lambda_1a = q_lambda(lam=lambda_1a, k=K_AIR)
        
        inlet = self.solve_inlet_geometry(q_lambda_1a, L_ku, a_1crit)
        u_mid = inlet.u_mid
        flow_coefficient = self.params.c_1a / u_mid

        stage_1 = calculate_section_diameters(
            D_ref=self.params.D_const,
            F=inlet.area,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )

        h_1 = (stage_1.tip - stage_1.hub) / 2

        alpha1_rad = math.radians(self.inlet_triangle.alpha_deg)
        self.inlet_triangle.beta_deg = math.degrees(math.atan(flow_coefficient / (1 - flow_coefficient * (1 / math.tan(alpha1_rad)))))


        self.inlet_triangle.w = math.sqrt(self.inlet_triangle.c_a ** 2 + (u_mid - self.inlet_triangle.c_u) ** 2)

        a_1 = local_speed_velocity(a_crit=a_1crit, k=K_AIR, lam=inlet.lambda_1)

        self.inlet_triangle.mach = self.inlet_triangle.w  / a_1

        par_L_c = inlet.L_ku_rel / flow_coefficient
        par_rho_c = self.params.reaction / flow_coefficient
        par_L_c_bt = 0.7 - 0.27 * par_rho_c + 0.16 * par_rho_c ** 2
        J = par_L_c / par_L_c_bt
        rotor_solidity = 0.225 + 0.275 * J + 0.5 * J ** 2


        z_rotor = self.params.h_rot_rel * rotor_solidity * math.pi * stage_1.mid / h_1
        z_rotor = math.ceil(z_rotor)
        h_rot_rel = h_1 * z_rotor / (rotor_solidity * math.pi * stage_1.mid)

        rotor_blade_chord = h_1 / h_rot_rel
        
        c_2u = u_mid * ((1 - self.params.reaction) + inlet.L_ku_rel / 2)
        c_2a = (self.params.c_1a + self.params.c_3a) / 2

        self.outlet_triangle = VelocityTriangle(c=math.sqrt(c_2a ** 2 + c_2u ** 2),
                                                c_a=c_2a,
                                                c_u=c_2u)
        a_2crit = velocity_critical(k=K_AIR, R=R_AIR, T=thermo.T_out)
        lambda_2 = self.outlet_triangle.c / a_2crit
        a_2 = local_speed_velocity(a_crit=a_2crit, k=K_AIR, lam=lambda_2)

        alpha2_rad = math.asin(self.outlet_triangle.c_a / self.outlet_triangle.c)
        self.outlet_triangle.mach = self.outlet_triangle.c / a_2
        self.outlet_triangle.alpha_deg = math.degrees(alpha2_rad)


        
        sigma_stator = thermo.p_out / thermo.p_2

        q_lambda_2 = q_lambda(lam=lambda_2, k=K_AIR)
        F_outlet = area_calc(G=self.params.G, T=thermo.T_out, s=S_AIR, p=thermo.p_2,
                            q=(q_lambda_2 * self.params.K_g * math.sin(alpha2_rad)))
        
        d_2_hub_rel = relative_diameter_hub(D_ref=self.params.D_const, F=F_outlet, mode_name=self.params.mode)

        stage_2 = calculate_section_diameters(
            D_ref=self.params.D_const,
            F=F_outlet,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )
        h_2 = (stage_2.tip - stage_2.hub) / 2
        self.outlet_triangle.w = math.sqrt(self.outlet_triangle.c_a ** 2 + (u_mid - self.outlet_triangle.c_u) ** 2)

        beta_2_rad = math.asin(self.outlet_triangle.c_a / self.outlet_triangle.w)
        self.outlet_triangle.beta_deg = math.degrees(beta_2_rad)


        # print(f'pi = {pi}')
        # print(f'u_mid = {u_mid}')
        # print(f'a_crit = {a_1crit}')
        # print(f'lambda 1a = {lambda_1a}')
        # print(f'q_lambda_a1 {q_lambda_1a}')

        # print(f'u_mid = {u_mid}')

        # print(f'коэфф расхода {flow_coefficient}')

        # print(f'Диаметры на входе {stage_1}')
        # print(f'высота лопатки РК {h_1}')
        # print(f'скорость звука на входе {a_1}')
        # print(f'Треугольник скоростей вход {self.inlet_triangle}')
        # print(f'безразмерный параметр J {J}')
        # print(f'густота решетки РК {rotor_solidity}')
        # print(f'число лопаток РК {z_rotor}')
        # print(f'относительная высота рабочей лопатки {h_rot_rel}')
        # print(f'хорда рабочих лопаток {rotor_blade_chord}')

        # print(f'lamda 2 = {lambda_2}')
        # print(f'скорость звука на выходе из РК {a_2}')
        # print(f'Треугольник выхода из РК {self.outlet_triangle}')
        # print(f'Площадь на выходе из РК {F_outlet}')
        # print(f'относительный диаметр втулки на выходе из РК {d_2_hub_rel}')
        # print(f'Диаметры на выходе из РК {stage_2}')
        # print(f'высота лопатки на выходе из РК {h_2}')

    def calculate_thermodynamics(self) -> StageThermodynamics:
        L_st = self.params.L_rel * self.params.u_tip ** 2
        T_out = self.params.T_in + L_st / (K_AIR / (K_AIR - 1) * R_AIR) # T_in - температура перед прошлой ступенью
        pi = (L_st * self.params.efficiency / (K_AIR / (K_AIR - 1) * R_AIR * self.params.T_in) + 1) ** (K_AIR / (K_AIR - 1))
        p_out = self.params.p_in * pi
        p_2 = self.params.p_in * (1 + L_st * self.params.eff_rotor / (K_AIR / (K_AIR - 1) * R_AIR * self.params.T_in)) ** (K_AIR / (K_AIR - 1))
        return StageThermodynamics(L_stage=L_st, T_in=self.params.T_in,
                                   T_out=T_out, p_in=self.params.p_in, p_2=p_2, p_out=p_out, pi=pi)
    
    def solve_inlet_geometry(self, q_lam_1a, L_ku, a_crit) -> InletSolution: 
        delta = 10
        alpha_deg = 90
        q_lambda_1 = q_lam_1a
        while delta > 1:
            alpha_old = alpha_deg
            alpha_rad = math.radians(alpha_old)
            F = area_calc(G=self.params.G, T=self.params.T_in, s=S_AIR, p=self.params.p_in,
                            q=(q_lambda_1 * self.params.K_g * math.sin(alpha_rad)))

            d_hub_rel = relative_diameter_hub(D_ref=self.params.D_const, F=F, mode_name=self.params.mode)

            u_mid = self.params.u_tip * math.sqrt((1 + d_hub_rel ** 2) / 2) # реализовать расчет при разных D=const

            L_ku_rel = L_ku / u_mid ** 2

            c_u = u_mid * ((1 - self.params.reaction) - L_ku_rel / 2)

            c = math.sqrt(self.params.c_1a ** 2 + c_u ** 2)

            lambda_1 = c / a_crit
            q_lambda_1 = q_lambda(lam=lambda_1, k=K_AIR)

            alpha_deg = math.degrees(math.asin(self.params.c_1a / c))
            delta = abs((alpha_deg - alpha_old) / alpha_deg) * 100
            
        self.inlet_triangle = VelocityTriangle(
            c=c,
            c_a=self.params.c_1a,
            c_u=c_u,
            alpha_deg=alpha_deg
        )
        return InletSolution(area=F, u_mid=u_mid, lambda_1=lambda_1, L_ku_rel=L_ku_rel)

if __name__ == '__main__':
    c1a_stage.append(180.8988)
    for i in range(1):
        params_stage = CompressorStageParameters(
            mode='tip',
            D_const= D_tip,

            alpha1_deg=alpha11_deg,
            L_rel=L_stage_rel[i],
            D_tip=D_tip,
            efficiency=eff_stage[i],
            eff_tip_rel=eff_rel[i],
            reaction=reaction_stage[i],
            c_1a=c1a_stage[i],
            c_3a = c1a_stage[i+1],
            h_rot_rel=h_rot_rel[i],
            G=G,
            d_hub_rel=d_hub_rel,
            lambda_in=lambda_in,
            u_tip=u_tip_1,
            T_in=T_in,
            p_in=p_in,
            D_1mid=D_1mid_1,
            D_1hub=D_1hub_1,
            K_g = 0.96,
            eff_rotor=0.93
        )
        com_stage = CompressorStage(params_stage)
        com_stage.calculate()
