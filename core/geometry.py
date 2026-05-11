import math
from core.geometry_models import SectionDiameters

def area_calc(G, T, s, p, q):
    return G * math.sqrt(T) / (s * p * q)

def diameter_from_area(D_ref, F, k):
    return math.sqrt(D_ref**2 + k * F / math.pi)

def calculate_section_diameters(D_ref, F, mode_name, MODES_D):
    mode = MODES_D[mode_name]

    return SectionDiameters(
        hub=diameter_from_area(D_ref, F, mode['hub']),
        mid=diameter_from_area(D_ref, F, mode['mid']),
        tip=diameter_from_area(D_ref, F, mode['tip']),
    )

def reconstruct_section_diameters(F, hub_to_tip_ratio):
    D_tip = math.sqrt(
        4 * F / (math.pi * (1 - hub_to_tip_ratio**2))
    )

    D_hub = hub_to_tip_ratio * D_tip

    D_mid = D_tip * math.sqrt(
        (1 + hub_to_tip_ratio**2) / 2
    )

    return SectionDiameters(
        hub=D_hub,
        mid=D_mid,
        tip=D_tip
    )

def relative_diameter_hub(D_ref, F, mode_name):
    K = {'tip': (4,0),
         'mid': (2,2),
         'hub': (0,4)}
    k1, k2 = K[mode_name]
    d = math.sqrt((D_ref ** 2 - k1 * F / math.pi) / (D_ref ** 2 + k2 * F / math.pi))
    return d
