import math
from core.diametrs_classes import SectionDiameters

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