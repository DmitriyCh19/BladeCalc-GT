from stages.hpt import HPTParameters, HighPressureTurbine


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

if __name__ == '__main__':
    hpt = HighPressureTurbine(hpt_params)

    hpt.calculate()

    print(hpt.geometry)
    print(hpt.n_hp)