from dataclasses import dataclass

from core.geometry_models import SectionDiameters, VelocityTriangle


@dataclass
class TurbineStageFlowParams:
    G_in: float
    G_out: float


@dataclass
class TurbineStageThermoParams:
    p_in: float
    pi: float
    T_in: float
    T_out: float
    L_stage: float


@dataclass
class TurbineStageKinematicParams:
    u_mid: float
    y: float
    reaction: float
    phi_cooling: float = 1.0
    psi_cooling: float = 1.0


@dataclass
class TurbineStageGeometryParams:
    mode: str
    inlet: SectionDiameters
    h_2_rotor: float


@dataclass
class BladeRowInput:
    axial_chord: float
    trailing_edge_radius_rel: float
    cooling_pitch_factor: float = 1.0


@dataclass
class TurbineStageBladeRowsParams:
    stator: BladeRowInput
    rotor: BladeRowInput


@dataclass
class TurbineStageParameters:
    flow: TurbineStageFlowParams
    thermo: TurbineStageThermoParams
    kinematics: TurbineStageKinematicParams
    geometry: TurbineStageGeometryParams
    blade_rows: TurbineStageBladeRowsParams


@dataclass
class TurbineStageWork:
    L_ad: float
    L_ad_stator: float
    L_ad_rotor: float
    L_stage: float
    c_ad: float
    lambda_ad: float
    p_2_static_ad: float
    y: float


@dataclass
class TurbineStatorOutlet:
    triangle: VelocityTriangle
    T_static: float
    p_static: float
    density: float
    lambda_1: float
    lambda_1t: float
    q_lambda_1: float
    area: float
    section: SectionDiameters
    h_stator_in: float
    h_stator_out: float
    K_reaction: float
    reaction_hub: float


@dataclass
class TurbineRotorOutlet:
    triangle: VelocityTriangle
    T_relative_in: float
    T_static: float
    p_relative_out: float
    p_static: float
    lambda_w2: float
    lambda_c2: float
    q_lambda_w2: float
    q_lambda_2: float
    area: float
    section: SectionDiameters


@dataclass
class BladeRowGridParams:
    t_opt_rel: float
    delta_t_opt_rel: float
    gamma_deg: float
    blade_chord: float
    trailing_edge_radius: float
    S_out_rel: float
    K_crit: float
    pitch: float
    pitch_with_cooling: float
    solidity: float
    blade_count: int


@dataclass
class TurbineStageBladeRows:
    stator_grid: BladeRowGridParams
    rotor_grid: BladeRowGridParams
    delta_alpha_deg: float
    delta_beta_deg: float
    K_stator: float
    K_rotor: float
    solidity_hub_rotor: float


@dataclass
class TurbineStageThermodynamics:
    L_ad: float
    L_ad_stator: float
    L_ad_rotor: float
    L_stage: float
    T_in: float
    T_1: float
    T_w1: float
    T_2: float
    T_out_target: float
    p_in: float
    p_1_static: float
    p_2_static: float
    p_w2: float
    p_out: float
    pi_stage: float


@dataclass
class TurbineStageGeometry:
    inlet: SectionDiameters
    stator_outlet: SectionDiameters
    rotor_outlet: SectionDiameters
    h_0_stator: float
    h_1_stator: float
    h_2_rotor: float
    F_1: float
    F_2: float


@dataclass
class TurbineStageLoading:
    reaction_mean: float
    reaction_hub: float
    K_reaction: float
    delta_alpha_deg: float
    delta_beta_deg: float
    K_stator: float
    K_rotor: float


@dataclass
class TurbineStageVelocity:
    stator_outlet: VelocityTriangle
    rotor_outlet: VelocityTriangle
    c_ad: float
    lambda_ad: float
    lambda_1: float
    lambda_1t: float
    lambda_w2: float
    lambda_c2: float


@dataclass
class TurbineStageResult:
    thermodynamics: TurbineStageThermodynamics
    geometry: TurbineStageGeometry
    velocity: TurbineStageVelocity
    loading: TurbineStageLoading
    stator_grid: BladeRowGridParams
    rotor_grid: BladeRowGridParams
    solidity_hub_rotor: float