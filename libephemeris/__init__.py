from .constants import *
from .exceptions import Error
from .time_utils import (
    swe_julday,
    swe_revjul,
    swe_deltat,
    swe_deltat_ex,
    date_conversion,
    day_of_week,
    utc_to_jd,
    jdet_to_utc,
    jdut1_to_utc,
    utc_time_zone,
    time_equ,
    lat_to_lmt,
    lmt_to_lat,
    sidtime,
    sidtime0,
)
from .planets import (
    swe_calc_ut,
    swe_calc,
    swe_calc_pctr,
    swe_get_ayanamsa_ut,
    swe_get_ayanamsa,
    swe_get_ayanamsa_ex,
    swe_get_ayanamsa_ex_ut,
    swe_get_ayanamsa_name,
    swe_set_sid_mode,
    swe_nod_aps,
    swe_nod_aps_ut,
    swe_get_orbital_elements,
    swe_get_orbital_elements_ut,
    swe_orbit_max_min_true_distance,
    swe_pheno,
    swe_pheno_ut,
    get_planet_name,
)
from .houses import (
    swe_houses,
    swe_houses_armc,
    swe_houses_armc_ex2,
    swe_houses_ex,
    swe_houses_ex2,
    swe_house_name,
    swe_house_pos,
    house_pos,
    gauquelin_sector,
    swe_gauquelin_sector,
)
from .state import (
    set_topo as swe_set_topo,
    set_ephe_path as swe_set_ephe_path,
    set_ephemeris_file as swe_set_ephemeris_file,
    set_jpl_file as swe_set_jpl_file,
    set_tid_acc as swe_set_tid_acc,
    get_tid_acc as swe_get_tid_acc,
    set_delta_t_userdef as swe_set_delta_t_userdef,
    get_delta_t_userdef as swe_get_delta_t_userdef,
    set_lapse_rate as swe_set_lapse_rate,
    get_lapse_rate as swe_get_lapse_rate,
    get_library_path as swe_get_library_path,
    get_current_file_data as swe_get_current_file_data,
    close as swe_close,
)
from .crossing import (
    swe_solcross_ut,
    swe_solcross,
    swe_mooncross_ut,
    swe_mooncross,
    swe_mooncross_node_ut,
    swe_mooncross_node,
    swe_cross_ut,
    swe_helio_cross_ut,
    swe_helio_cross,
)
from .eclipse import (
    sol_eclipse_when_glob,
    swe_sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    swe_sol_eclipse_when_loc,
    sol_eclipse_where,
    swe_sol_eclipse_where,
    sol_eclipse_how,
    swe_sol_eclipse_how,
    lun_eclipse_when,
    swe_lun_eclipse_when,
    lun_eclipse_when_loc,
    swe_lun_eclipse_when_loc,
    lun_eclipse_how,
    swe_lun_eclipse_how,
    lun_occult_when_glob,
    swe_lun_occult_when_glob,
    lun_occult_when_loc,
    swe_lun_occult_when_loc,
    lun_occult_where,
    swe_lun_occult_where,
    rise_trans,
    swe_rise_trans,
    rise_trans_true_hor,
    swe_rise_trans_true_hor,
    heliacal_ut,
    swe_heliacal_ut,
    heliacal_pheno_ut,
    swe_heliacal_pheno_ut,
    vis_limit_mag,
    swe_vis_limit_mag,
)
from .utils import (
    degnorm,
    radnorm,
    difdeg2n,
    difdegn,
    difrad2n,
    difcs2n,
    difcsn,
    csnorm,
    csroundsec,
    cs2degstr,
    cs2lonlatstr,
    cs2timestr,
    deg_midp,
    rad_midp,
    d2l,
    split_deg,
    swe_calc_angles,
    cotrans,
    cotrans_sp,
    azalt,
    azalt_rev,
    refrac,
    refrac_extended,
    SE_ECL2HOR,
    SE_EQU2HOR,
    SE_HOR2ECL,
    SE_HOR2EQU,
    SE_TRUE_TO_APP,
    SE_APP_TO_TRUE,
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ZODIACAL,
    SPLIT_DEG_NAKSHATRA,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_KEEP_DEG,
)
from .fixed_stars import (
    swe_fixstar_ut,
    swe_fixstar,
    swe_fixstar2,
    swe_fixstar2_ut,
    swe_fixstar_mag,
    swe_fixstar2_mag,
)
from .context import EphemerisContext  # NEW: Thread-safe context API


# =============================================================================
# PYSWISSEPH-COMPATIBLE FUNCTION ALIASES (without swe_ prefix)
# =============================================================================
# pyswisseph uses function names without the swe_ prefix
# These aliases provide 100% API compatibility with pyswisseph

# Time functions
julday = swe_julday
revjul = swe_revjul
deltat = swe_deltat
deltat_ex = swe_deltat_ex
# date_conversion already uses snake_case, no alias needed

# Planet calculation
calc_ut = swe_calc_ut
calc = swe_calc
calc_pctr = swe_calc_pctr
nod_aps = swe_nod_aps
nod_aps_ut = swe_nod_aps_ut
get_orbital_elements = swe_get_orbital_elements
get_orbital_elements_ut = swe_get_orbital_elements_ut
orbit_max_min_true_distance = swe_orbit_max_min_true_distance
pheno = swe_pheno
pheno_ut = swe_pheno_ut

# Houses
houses = swe_houses
houses_armc = swe_houses_armc
houses_armc_ex2 = swe_houses_armc_ex2
houses_ex = swe_houses_ex
houses_ex2 = swe_houses_ex2
house_name = swe_house_name
# house_pos is already the main function name (matching pyswisseph)

# Ayanamsa (sidereal)
get_ayanamsa_ut = swe_get_ayanamsa_ut
get_ayanamsa = swe_get_ayanamsa
get_ayanamsa_ex = swe_get_ayanamsa_ex
get_ayanamsa_ex_ut = swe_get_ayanamsa_ex_ut
get_ayanamsa_name = swe_get_ayanamsa_name
set_sid_mode = swe_set_sid_mode

# Observer location
set_topo = swe_set_topo
set_ephe_path = swe_set_ephe_path
set_ephemeris_file = swe_set_ephemeris_file
set_jpl_file = swe_set_jpl_file

# Tidal acceleration for Delta T
set_tid_acc = swe_set_tid_acc
get_tid_acc = swe_get_tid_acc

# User-defined Delta T
set_delta_t_userdef = swe_set_delta_t_userdef
get_delta_t_userdef = swe_get_delta_t_userdef

# Lapse rate for refraction calculations
set_lapse_rate = swe_set_lapse_rate
get_lapse_rate = swe_get_lapse_rate

# Library path
get_library_path = swe_get_library_path

# Current file data
get_current_file_data = swe_get_current_file_data

# Close and cleanup
close = swe_close

# Fixed Stars
fixstar_ut = swe_fixstar_ut
fixstar = swe_fixstar
fixstar2 = swe_fixstar2
fixstar2_ut = swe_fixstar2_ut
fixstar_mag = swe_fixstar_mag
fixstar2_mag = swe_fixstar2_mag

# Crossings
solcross_ut = swe_solcross_ut
solcross = swe_solcross
mooncross_ut = swe_mooncross_ut
mooncross = swe_mooncross
mooncross_node_ut = swe_mooncross_node_ut
mooncross_node = swe_mooncross_node
helio_cross_ut = swe_helio_cross_ut
helio_cross = swe_helio_cross


# Helper for Arabic parts
from .arabic_parts import calc_all_arabic_parts

# Constants (planet IDs, flags, sidereal modes)
from .constants import *

__version__ = "0.1.8"
__author__ = "Giacomo Battaglia"
__license__ = "LGPL-3.0"

__all__ = [
    # Exceptions
    "Error",
    # Thread-safe Context API
    "EphemerisContext",
    # Time functions (both swe_ and non-prefixed aliases)
    "swe_julday",
    "julday",
    "swe_revjul",
    "revjul",
    "swe_deltat",
    "deltat",
    "swe_deltat_ex",
    "deltat_ex",
    "date_conversion",
    "day_of_week",
    "utc_to_jd",
    "jdet_to_utc",
    "jdut1_to_utc",
    "utc_time_zone",
    "time_equ",
    "lat_to_lmt",
    "lmt_to_lat",
    "sidtime",
    "sidtime0",
    # Planet calculation
    "swe_calc_ut",
    "calc_ut",
    "swe_calc",
    "calc",
    "swe_calc_pctr",
    "calc_pctr",
    "swe_nod_aps",
    "nod_aps",
    "swe_nod_aps_ut",
    "nod_aps_ut",
    "swe_get_orbital_elements",
    "get_orbital_elements",
    "swe_get_orbital_elements_ut",
    "get_orbital_elements_ut",
    "swe_orbit_max_min_true_distance",
    "orbit_max_min_true_distance",
    # Planetary phenomena
    "swe_pheno",
    "pheno",
    "swe_pheno_ut",
    "pheno_ut",
    # Houses
    "swe_houses",
    "houses",
    "swe_houses_armc",
    "houses_armc",
    "swe_houses_armc_ex2",
    "houses_armc_ex2",
    "swe_houses_ex",
    "houses_ex",
    "swe_houses_ex2",
    "houses_ex2",
    "swe_house_name",
    "house_name",
    "swe_house_pos",
    "house_pos",
    "gauquelin_sector",
    "swe_gauquelin_sector",
    # Ayanamsa (sidereal)
    "swe_set_sid_mode",
    "set_sid_mode",
    "swe_get_ayanamsa_ut",
    "get_ayanamsa_ut",
    "swe_get_ayanamsa",
    "get_ayanamsa",
    "swe_get_ayanamsa_ex",
    "get_ayanamsa_ex",
    "swe_get_ayanamsa_ex_ut",
    "get_ayanamsa_ex_ut",
    "swe_get_ayanamsa_name",
    "get_ayanamsa_name",
    # Observer location
    "swe_set_topo",
    "set_topo",
    "swe_set_ephe_path",
    "set_ephe_path",
    "swe_set_ephemeris_file",
    "set_ephemeris_file",
    "swe_set_jpl_file",
    "set_jpl_file",
    # Tidal acceleration
    "swe_set_tid_acc",
    "set_tid_acc",
    "swe_get_tid_acc",
    "get_tid_acc",
    # User-defined Delta T
    "swe_set_delta_t_userdef",
    "set_delta_t_userdef",
    "swe_get_delta_t_userdef",
    "get_delta_t_userdef",
    # Lapse rate for refraction
    "swe_set_lapse_rate",
    "set_lapse_rate",
    "swe_get_lapse_rate",
    "get_lapse_rate",
    # Library path
    "swe_get_library_path",
    "get_library_path",
    # Current file data
    "swe_get_current_file_data",
    "get_current_file_data",
    # Close and cleanup
    "swe_close",
    "close",
    # Crossings
    "swe_solcross_ut",
    "solcross_ut",
    "swe_solcross",
    "solcross",
    "swe_mooncross_ut",
    "mooncross_ut",
    "swe_mooncross",
    "mooncross",
    "swe_mooncross_node_ut",
    "mooncross_node_ut",
    "swe_mooncross_node",
    "mooncross_node",
    "swe_cross_ut",
    "swe_helio_cross_ut",
    "helio_cross_ut",
    "swe_helio_cross",
    "helio_cross",
    # Eclipses
    "sol_eclipse_when_glob",
    "swe_sol_eclipse_when_glob",
    "sol_eclipse_when_loc",
    "swe_sol_eclipse_when_loc",
    "sol_eclipse_where",
    "swe_sol_eclipse_where",
    "sol_eclipse_how",
    "swe_sol_eclipse_how",
    "lun_eclipse_when",
    "swe_lun_eclipse_when",
    "lun_eclipse_when_loc",
    "swe_lun_eclipse_when_loc",
    "lun_eclipse_how",
    "swe_lun_eclipse_how",
    "lun_occult_when_glob",
    "swe_lun_occult_when_glob",
    "lun_occult_when_loc",
    "swe_lun_occult_when_loc",
    "lun_occult_where",
    "swe_lun_occult_where",
    # Rise/Set/Transit
    "rise_trans",
    "swe_rise_trans",
    "rise_trans_true_hor",
    "swe_rise_trans_true_hor",
    # Heliacal events
    "heliacal_ut",
    "swe_heliacal_ut",
    "heliacal_pheno_ut",
    "swe_heliacal_pheno_ut",
    # Utilities
    "degnorm",
    "radnorm",
    "difdeg2n",
    "difdegn",
    "difrad2n",
    "difcs2n",
    "difcsn",
    "csnorm",
    "csroundsec",
    "cs2degstr",
    "cs2lonlatstr",
    "cs2timestr",
    "deg_midp",
    "rad_midp",
    "d2l",
    "split_deg",
    "swe_calc_angles",
    "cotrans",
    "cotrans_sp",
    "azalt",
    "azalt_rev",
    "refrac",
    "refrac_extended",
    "SE_ECL2HOR",
    "SE_EQU2HOR",
    "SE_HOR2ECL",
    "SE_HOR2EQU",
    "SE_TRUE_TO_APP",
    "SE_APP_TO_TRUE",
    "SPLIT_DEG_ROUND_SEC",
    "SPLIT_DEG_ROUND_MIN",
    "SPLIT_DEG_ROUND_DEG",
    "SPLIT_DEG_ZODIACAL",
    "SPLIT_DEG_NAKSHATRA",
    "SPLIT_DEG_KEEP_SIGN",
    "SPLIT_DEG_KEEP_DEG",
    # Helpers
    "calc_all_arabic_parts",
    "get_planet_name",
    # Fixed Stars
    "swe_fixstar_ut",
    "fixstar_ut",
    "swe_fixstar",
    "fixstar",
    "swe_fixstar2",
    "fixstar2",
    "swe_fixstar2_ut",
    "fixstar2_ut",
    "swe_fixstar_mag",
    "fixstar_mag",
    "swe_fixstar2_mag",
    "fixstar2_mag",
]
