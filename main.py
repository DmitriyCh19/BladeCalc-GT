from machine.hpt import HPTParameters, HighPressureTurbine
from machine.hpc import HPCParameters, HighPressureCompressor
from machine.compressor_stage import CompressorStageParameters, CompressorStage


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

hpt_params = HPTParameters(
    mode='mid',

    T_gas=1665,
    p_gas=2199900,

    T_out=1327.36,
    p_out=765440,

    T_cooling_air=796.99,

    z=1,
    efficiency=0.88,
    y=0.445,
    theta=0.58,

    G_in=66.496,
    G_cooling_rel=0.06,
    G_hpc_in=72.435,

    density=8100,
    blade_force_coeff=0.7,
    sigma_allow=380,
    k_sigma=1.8,

    alpha_2_deg=65,

    lambda_2=0.497,
    lambda_gas=0.25,

    L_stage=393673,
)

hpc_params = HPCParameters(
    mode='tip',

    T_in=437.97,
    T_out=796.99,
    p_out=2.3402e6,
    pi=6.58,
    d_hub_in_rel=0.74,
    K_g=0.96,
    efficiency=0.87,
    lambda_in=0.58,
    lambda_out=0.35,

    G_in=72.435,
    G_cooling_rel=0.1
)

if __name__ == '__main__':
    K_gg = 0.4

    hpt = HighPressureTurbine(hpt_params)
    hpc = HighPressureCompressor(hpc_params)

    hpt.calculate()
    hpc.calculate(n=hpt.n_hp)
    z_hpc = (hpt.D_mean / hpc.D_mean) ** 2 * hpt.params.z / K_gg ** 2

    print(hpt.geometry)
    print(f'Относительная высота лопатки = {hpt.blade_height_rel}')
    print(f'Высота лопатки на выходе из ТВД ≥0.025м = {hpt.h_blade_out}')
    print(hpc.geometry)
    print(f'c_out КВД 150...180 = {hpc.c_out}')
    print(f'u внешний диаметр КВД 1 ступень ≤ 550 м/с = {hpc.u_tip_in} ')
    print(z_hpc)

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
    
    com_stage.print_stage_result()