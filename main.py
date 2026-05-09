from machine.hpt import HPTParameters, HighPressureTurbine
from machine.hpc import HPCParameters, HighPressureCompressor

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