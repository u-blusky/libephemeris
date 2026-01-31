from .constants import *
from .exceptions import (
    # Base error class (pyswisseph compatible)
    Error,
    # Category: Input validation errors
    InputValidationError,
    CoordinateError,
    InvalidBodyError,
    # Category: Data not found errors
    DataNotFoundError,
    UnknownBodyError,
    StarNotFoundError,
    SPKNotFoundError,
    # Category: Calculation errors
    CalculationError,
    PolarCircleError,
    EphemerisRangeError,
    ConvergenceError,
    # Category: Configuration errors
    ConfigurationError,
    # Validation helpers
    validate_latitude,
    validate_longitude,
    validate_coordinates,
)
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
    # TAI (International Atomic Time) functions
    TT_TAI_OFFSET_SECONDS,
    TT_TAI_OFFSET_DAYS,
    get_tai_utc_for_jd,
    utc_to_tai_jd,
    tai_jd_to_utc,
    tt_to_tai_jd,
    tai_to_tt_jd,
    tai_to_utc_jd,
    utc_jd_to_tai,
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
    # Elongation helper functions
    get_elongation_from_sun,
    get_signed_elongation,
    is_morning_star,
    is_evening_star,
    get_elongation_type,
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
    # Polar latitude handling
    swe_houses_with_fallback,
    swe_houses_armc_with_fallback,
    get_polar_latitude_threshold,
    # Extreme latitude handling
    get_extreme_latitude_info,
    EXTREME_LATITUDE_THRESHOLD,
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
    set_auto_spk_download,
    get_auto_spk_download,
    set_spk_cache_dir,
    get_spk_cache_dir,
    set_spk_date_padding,
    get_spk_date_padding,
    # IERS Delta T configuration
    set_iers_delta_t_enabled,
    get_iers_delta_t_enabled,
)
from .iers_data import (
    # IERS data download functions
    download_iers_finals,
    download_leap_seconds,
    download_delta_t_data,
    load_iers_data,
    # IERS Delta T lookup functions
    get_observed_delta_t,
    get_observed_delta_t_data_range,
    is_observed_delta_t_available,
    get_delta_t_iers,
    get_ut1_utc,
    get_tai_utc,
    get_iers_data_range,
    is_iers_data_available,
    # IERS cache management
    clear_iers_cache,
    delete_iers_cache_files,
    get_iers_cache_info,
    # IERS configuration
    set_iers_cache_dir,
    get_iers_cache_dir,
    set_iers_auto_download,
    get_iers_auto_download,
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
    BesselianElements,
    sol_eclipse_max_time,
    swe_sol_eclipse_max_time,
    sol_eclipse_when_glob,
    swe_sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    swe_sol_eclipse_when_loc,
    sol_eclipse_where,
    swe_sol_eclipse_where,
    sol_eclipse_how,
    swe_sol_eclipse_how,
    sol_eclipse_how_details,
    swe_sol_eclipse_how_details,
    # Eclipse magnitude at location
    sol_eclipse_magnitude_at_loc,
    swe_sol_eclipse_magnitude_at_loc,
    # Eclipse obscuration at location
    sol_eclipse_obscuration_at_loc,
    swe_sol_eclipse_obscuration_at_loc,
    lun_eclipse_when,
    swe_lun_eclipse_when,
    lun_eclipse_when_loc,
    swe_lun_eclipse_when_loc,
    lun_eclipse_how,
    swe_lun_eclipse_how,
    # Lunar eclipse umbral magnitude
    lun_eclipse_umbral_magnitude,
    swe_lun_eclipse_umbral_magnitude,
    # Lunar eclipse penumbral magnitude
    lun_eclipse_penumbral_magnitude,
    swe_lun_eclipse_penumbral_magnitude,
    # Lunar eclipse gamma parameter
    lun_eclipse_gamma,
    swe_lun_eclipse_gamma,
    lun_occult_when_glob,
    swe_lun_occult_when_glob,
    lun_occult_when_loc,
    swe_lun_occult_when_loc,
    lun_occult_where,
    swe_lun_occult_where,
    # Planetary occultations
    planet_occult_when_glob,
    swe_planet_occult_when_glob,
    planet_occult_when_loc,
    swe_planet_occult_when_loc,
    rise_trans,
    swe_rise_trans,
    rise_trans_true_hor,
    swe_rise_trans_true_hor,
    calc_besselian_x,
    calc_besselian_y,
    calc_besselian_d,
    calc_besselian_l1,
    calc_besselian_l2,
    calc_besselian_mu,
    calc_besselian_dx_dt,
    calc_besselian_dy_dt,
    calc_besselian_dd_dt,
    calc_besselian_dl1_dt,
    calc_besselian_dl2_dt,
    calc_besselian_dmu_dt,
    interpolate_besselian_elements,
    calc_eclipse_first_contact_c1,
    calc_eclipse_second_contact_c2,
    calc_eclipse_third_contact_c3,
    calc_eclipse_fourth_contact_c4,
    calc_lunar_eclipse_penumbral_first_contact_p1,
    calc_lunar_eclipse_penumbral_fourth_contact_p4,
    calc_lunar_eclipse_umbral_first_contact_u1,
    calc_lunar_eclipse_umbral_second_contact_u2,
    calc_lunar_eclipse_umbral_third_contact_u3,
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    # Eclipse duration functions
    calc_solar_eclipse_duration,
    calc_lunar_eclipse_total_duration,
    calc_lunar_eclipse_umbral_duration,
    # Eclipse path width calculation
    calc_eclipse_path_width,
    swe_calc_eclipse_path_width,
    # Eclipse central line coordinates
    calc_eclipse_central_line,
    swe_calc_eclipse_central_line,
    # Eclipse northern limit coordinates
    calc_eclipse_northern_limit,
    swe_calc_eclipse_northern_limit,
    # Eclipse southern limit coordinates
    calc_eclipse_southern_limit,
    swe_calc_eclipse_southern_limit,
    # Saros series calculation
    get_saros_number,
    SAROS_CYCLE_DAYS,
    # Inex series calculation
    get_inex_number,
    INEX_CYCLE_DAYS,
)
from .heliacal import (
    heliacal_ut,
    swe_heliacal_ut,
    heliacal_pheno_ut,
    swe_heliacal_pheno_ut,
    vis_limit_mag,
    swe_vis_limit_mag,
    is_inner_planet,
    is_fixed_star,
    INNER_PLANETS,
)
from .extinction import (
    # Extinction calculation functions
    calc_airmass,
    calc_extinction_coefficient,
    calc_extinction_magnitude,
    calc_simple_extinction,
    apparent_magnitude_with_extinction,
    get_extinction_for_heliacal,
    # Individual component functions
    calc_rayleigh_coefficient,
    calc_aerosol_coefficient,
    calc_ozone_coefficient,
    calc_water_vapor_coefficient,
    # Data class for extinction components
    ExtinctionCoefficients,
    # Constants
    WAVELENGTH_U,
    WAVELENGTH_B,
    WAVELENGTH_V,
    WAVELENGTH_R,
    WAVELENGTH_I,
    # Twilight sky brightness model
    TwilightSkyBrightness,
    get_twilight_phase,
    calc_twilight_sky_brightness,
    calc_twilight_brightness_simple,
    calc_limiting_magnitude_twilight,
    # Twilight phase constants
    TWILIGHT_CIVIL_START,
    TWILIGHT_CIVIL_END,
    TWILIGHT_NAUTICAL_END,
    TWILIGHT_ASTRONOMICAL_END,
    DARK_SKY_BRIGHTNESS_V,
    # Schaefer visibility threshold model
    VisibilityResult,
    calc_eye_adaptation_state,
    calc_contrast_threshold,
    calc_visibility_threshold,
    is_object_visible,
    calc_limiting_magnitude_for_sky,
    OBSERVER_SKILL_INEXPERIENCED,
    OBSERVER_SKILL_AVERAGE,
    OBSERVER_SKILL_EXPERIENCED,
    OBSERVER_SKILL_EXPERT,
    EXPERIENCE_FACTORS,
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
    propagate_proper_motion,
)
from .context import EphemerisContext  # NEW: Thread-safe context API
from .spk import (  # SPK kernel support for high-precision minor body calculations
    download_spk,
    register_spk_body,
    unregister_spk_body,
    get_spk_body_info,
    list_spk_bodies,
    get_spk_coverage,
    download_and_register_spk,
)
from . import spk_auto  # Automatic SPK download and caching
from .minor_bodies import (  # Generic asteroid lookup by number
    calc_asteroid_by_number,
    fetch_orbital_elements_from_sbdb,
    clear_asteroid_elements_cache,
    get_asteroid_number,
    clear_asteroid_name_cache,
)
from .hypothetical import (  # Hamburg School Uranian planets
    calc_cupido,
    calc_hades,
    calc_zeus,
    calc_kronos,
    calc_apollon,
    calc_admetos,
    calc_vulkanus,
    calc_poseidon,
    calc_transpluto,  # Transpluto (Isis) - hypothetical trans-Plutonian planet
    calc_vulcan,  # Vulcan - hypothetical intramercurial planet
    calc_waldemath,  # Waldemath Moon - hypothetical second moon of Earth (seorbel.txt #18)
    calc_proserpina,  # Proserpina - hypothetical trans-Plutonian planet
    calc_planet_x_pickering,  # Pickering's Planet O/X prediction (1919)
    calc_white_moon_position,  # White Moon (Selena) - opposite to Black Moon Lilith
    calc_uranian_planet,
    calc_hypothetical_position,
    URANIAN_KEPLERIAN_ELEMENTS,
    TRANSPLUTO_KEPLERIAN_ELEMENTS,
    VULCAN_ELEMENTS,  # Vulcan orbital elements from seorbel.txt
    WALDEMATH_ELEMENTS,  # Waldemath Moon orbital elements from seorbel.txt #18
    PICKERING_PLANET_X_ELEMENTS,  # Pickering's Planet O/X orbital elements
    # seorbel.txt parser for custom hypothetical bodies
    parse_seorbel,
    get_bundled_seorbel_path,
    load_bundled_seorbel,
    SeorbelElements,
    TPolynomial,
    get_seorbel_body_by_name,
    calc_seorbel_position,
    SE_CUPIDO as SE_CUPIDO_HYPO,  # Alias to avoid conflict with constants.py
    SE_HADES as SE_HADES_HYPO,  # Alias to avoid conflict with constants.py
    SE_ZEUS as SE_ZEUS_HYPO,  # Alias to avoid conflict with constants.py
    SE_KRONOS as SE_KRONOS_HYPO,  # Alias to avoid conflict with constants.py
    SE_APOLLON as SE_APOLLON_HYPO,  # Alias to avoid conflict with constants.py
    SE_ADMETOS as SE_ADMETOS_HYPO,  # Alias to avoid conflict with constants.py
    SE_VULKANUS as SE_VULKANUS_HYPO,  # Alias to avoid conflict with constants.py
    SE_POSEIDON as SE_POSEIDON_HYPO,  # Alias to avoid conflict with constants.py
    SE_ISIS as SE_ISIS_HYPO,  # Alias to avoid conflict with constants.py
    SE_TRANSPLUTO as SE_TRANSPLUTO_HYPO,  # Alias to avoid conflict with constants.py
    SE_VULCAN as SE_VULCAN_HYPO,  # Alias to avoid conflict with constants.py
    SE_WALDEMATH as SE_WALDEMATH_HYPO,  # Alias to avoid conflict with constants.py
    SE_WHITE_MOON as SE_WHITE_MOON_HYPO,  # Alias to avoid conflict with constants.py
    SE_PROSERPINA as SE_PROSERPINA_HYPO,  # Alias to avoid conflict with constants.py
    SE_PLANET_X_PICKERING as SE_PLANET_X_PICKERING_HYPO,  # Alias to avoid conflict with constants.py
)

# REBOUND/ASSIST n-body integration (optional dependency)
# Lazy import to avoid requiring rebound/assist at module load time
from .rebound_integration import (
    ReboundIntegrator,
    AssistEphemConfig,
    PropagationResult,
    check_rebound_available,
    check_assist_available,
    get_rebound_version,
    get_assist_version,
    propagate_orbit_rebound,
    propagate_orbit_assist,
    propagate_trajectory,
    compare_with_keplerian,
)

from .planetary_moons import (  # Planetary moons (Galilean moons, Titan, etc.)
    register_moon_spk,
    unregister_moon_spk,
    list_registered_moons,
    get_moon_name,
    is_planetary_moon,
    get_moon_coverage,
    calc_moon_position,
    close_moon_kernels,
    # NAIF IDs for planetary moons
    NAIF_IO,
    NAIF_EUROPA,
    NAIF_GANYMEDE,
    NAIF_CALLISTO,
    NAIF_MIMAS,
    NAIF_ENCELADUS,
    NAIF_TETHYS,
    NAIF_DIONE,
    NAIF_RHEA,
    NAIF_TITAN,
    NAIF_HYPERION,
    NAIF_IAPETUS,
    NAIF_MIRANDA,
    NAIF_ARIEL,
    NAIF_UMBRIEL,
    NAIF_TITANIA,
    NAIF_OBERON,
    NAIF_TRITON,
    NAIF_PHOBOS,
    NAIF_DEIMOS,
    NAIF_CHARON,
    # Moon ID to NAIF mapping
    MOON_NAIF_MAP,
    MOON_NAMES,
    MOON_PARENT_MAP,
)

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

__version__ = "0.2.0"
__author__ = "Giacomo Battaglia"
__license__ = "LGPL-3.0"

__all__ = [
    # Exceptions - Base
    "Error",
    # Exceptions - Input Validation Category
    "InputValidationError",
    "CoordinateError",
    "InvalidBodyError",
    # Exceptions - Data Not Found Category
    "DataNotFoundError",
    "UnknownBodyError",
    "StarNotFoundError",
    "SPKNotFoundError",
    # Exceptions - Calculation Category
    "CalculationError",
    "PolarCircleError",
    "EphemerisRangeError",
    "ConvergenceError",
    # Exceptions - Configuration Category
    "ConfigurationError",
    # Coordinate validation
    "validate_latitude",
    "validate_longitude",
    "validate_coordinates",
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
    # TAI (International Atomic Time) functions
    "TT_TAI_OFFSET_SECONDS",
    "TT_TAI_OFFSET_DAYS",
    "get_tai_utc_for_jd",
    "utc_to_tai_jd",
    "tai_jd_to_utc",
    "tt_to_tai_jd",
    "tai_to_tt_jd",
    "tai_to_utc_jd",
    "utc_jd_to_tai",
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
    # Elongation helper functions
    "get_elongation_from_sun",
    "get_signed_elongation",
    "is_morning_star",
    "is_evening_star",
    "get_elongation_type",
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
    "sol_eclipse_max_time",
    "swe_sol_eclipse_max_time",
    "sol_eclipse_when_glob",
    "swe_sol_eclipse_when_glob",
    "sol_eclipse_when_loc",
    "swe_sol_eclipse_when_loc",
    "sol_eclipse_where",
    "swe_sol_eclipse_where",
    "sol_eclipse_how",
    "swe_sol_eclipse_how",
    "sol_eclipse_how_details",
    "swe_sol_eclipse_how_details",
    # Eclipse magnitude at location
    "sol_eclipse_magnitude_at_loc",
    "swe_sol_eclipse_magnitude_at_loc",
    # Eclipse obscuration at location
    "sol_eclipse_obscuration_at_loc",
    "swe_sol_eclipse_obscuration_at_loc",
    "lun_eclipse_when",
    "swe_lun_eclipse_when",
    "lun_eclipse_when_loc",
    "swe_lun_eclipse_when_loc",
    "lun_eclipse_how",
    "swe_lun_eclipse_how",
    # Lunar eclipse umbral magnitude
    "lun_eclipse_umbral_magnitude",
    "swe_lun_eclipse_umbral_magnitude",
    # Lunar eclipse penumbral magnitude
    "lun_eclipse_penumbral_magnitude",
    "swe_lun_eclipse_penumbral_magnitude",
    "lun_occult_when_glob",
    "swe_lun_occult_when_glob",
    "lun_occult_when_loc",
    "swe_lun_occult_when_loc",
    "lun_occult_where",
    "swe_lun_occult_where",
    # Besselian elements
    "calc_besselian_x",
    "calc_besselian_y",
    "calc_besselian_d",
    "calc_besselian_l1",
    "calc_besselian_l2",
    "calc_besselian_mu",
    # Besselian element time derivatives
    "calc_besselian_dx_dt",
    "calc_besselian_dy_dt",
    "calc_besselian_dd_dt",
    "calc_besselian_dl1_dt",
    "calc_besselian_dl2_dt",
    "calc_besselian_dmu_dt",
    # Solar eclipse contact points
    "calc_eclipse_first_contact_c1",
    "calc_eclipse_second_contact_c2",
    "calc_eclipse_third_contact_c3",
    "calc_eclipse_fourth_contact_c4",
    # Lunar eclipse penumbral contact points
    "calc_lunar_eclipse_penumbral_first_contact_p1",
    "calc_lunar_eclipse_penumbral_fourth_contact_p4",
    # Lunar eclipse umbral contact points
    "calc_lunar_eclipse_umbral_first_contact_u1",
    "calc_lunar_eclipse_umbral_second_contact_u2",
    "calc_lunar_eclipse_umbral_third_contact_u3",
    "calc_lunar_eclipse_umbral_fourth_contact_u4",
    # Eclipse duration functions
    "calc_solar_eclipse_duration",
    "calc_lunar_eclipse_total_duration",
    "calc_lunar_eclipse_umbral_duration",
    # Eclipse path width calculation
    "calc_eclipse_path_width",
    "swe_calc_eclipse_path_width",
    # Eclipse central line coordinates
    "calc_eclipse_central_line",
    "swe_calc_eclipse_central_line",
    # Saros series calculation
    "get_saros_number",
    "SAROS_CYCLE_DAYS",
    # Inex series calculation
    "get_inex_number",
    "INEX_CYCLE_DAYS",
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
    "vis_limit_mag",
    "swe_vis_limit_mag",
    "is_inner_planet",
    "is_fixed_star",
    "INNER_PLANETS",
    # Atmospheric extinction
    "calc_airmass",
    "calc_extinction_coefficient",
    "calc_extinction_magnitude",
    "calc_simple_extinction",
    "apparent_magnitude_with_extinction",
    "get_extinction_for_heliacal",
    "calc_rayleigh_coefficient",
    "calc_aerosol_coefficient",
    "calc_ozone_coefficient",
    "calc_water_vapor_coefficient",
    "ExtinctionCoefficients",
    "WAVELENGTH_U",
    "WAVELENGTH_B",
    "WAVELENGTH_V",
    "WAVELENGTH_R",
    "WAVELENGTH_I",
    # Twilight sky brightness model
    "TwilightSkyBrightness",
    "get_twilight_phase",
    "calc_twilight_sky_brightness",
    "calc_twilight_brightness_simple",
    "calc_limiting_magnitude_twilight",
    "TWILIGHT_CIVIL_START",
    "TWILIGHT_CIVIL_END",
    "TWILIGHT_NAUTICAL_END",
    "TWILIGHT_ASTRONOMICAL_END",
    "DARK_SKY_BRIGHTNESS_V",
    # Schaefer visibility threshold model
    "VisibilityResult",
    "calc_eye_adaptation_state",
    "calc_contrast_threshold",
    "calc_visibility_threshold",
    "is_object_visible",
    "calc_limiting_magnitude_for_sky",
    "OBSERVER_SKILL_INEXPERIENCED",
    "OBSERVER_SKILL_AVERAGE",
    "OBSERVER_SKILL_EXPERIENCED",
    "OBSERVER_SKILL_EXPERT",
    "EXPERIENCE_FACTORS",
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
    # Proper motion propagation
    "propagate_proper_motion",
    # SPK kernel support
    "download_spk",
    "register_spk_body",
    "unregister_spk_body",
    "get_spk_body_info",
    "list_spk_bodies",
    "get_spk_coverage",
    "download_and_register_spk",
    # Automatic SPK download module
    "spk_auto",
    # Auto SPK download configuration
    "set_auto_spk_download",
    "get_auto_spk_download",
    # SPK cache directory configuration
    "set_spk_cache_dir",
    "get_spk_cache_dir",
    # SPK date padding configuration
    "set_spk_date_padding",
    "get_spk_date_padding",
    # SPK body name mapping
    "SPK_BODY_NAME_MAP",
    "get_horizons_id",
    "get_naif_id_from_ipl",
    "get_spk_body_info_from_map",
    # IERS Delta T configuration
    "set_iers_delta_t_enabled",
    "get_iers_delta_t_enabled",
    # IERS data download functions
    "download_iers_finals",
    "download_leap_seconds",
    "download_delta_t_data",
    "load_iers_data",
    # IERS Delta T lookup functions
    "get_observed_delta_t",
    "get_observed_delta_t_data_range",
    "is_observed_delta_t_available",
    "get_delta_t_iers",
    "get_ut1_utc",
    "get_tai_utc",
    "get_iers_data_range",
    "is_iers_data_available",
    # IERS cache management
    "clear_iers_cache",
    "delete_iers_cache_files",
    "get_iers_cache_info",
    # IERS cache directory configuration
    "set_iers_cache_dir",
    "get_iers_cache_dir",
    # IERS auto-download configuration
    "set_iers_auto_download",
    "get_iers_auto_download",
    # Asteroid name lookup
    "get_asteroid_number",
    "clear_asteroid_name_cache",
    "calc_asteroid_by_number",
    "fetch_orbital_elements_from_sbdb",
    "clear_asteroid_elements_cache",
    # Hypothetical bodies (Hamburg School Uranian planets)
    "calc_cupido",
    "calc_hades",
    "calc_zeus",
    "calc_kronos",
    "calc_apollon",
    "calc_admetos",
    "calc_vulkanus",
    "calc_poseidon",
    "calc_transpluto",  # Transpluto (Isis) - hypothetical trans-Plutonian planet
    "calc_vulcan",  # Vulcan - hypothetical intramercurial planet
    "calc_waldemath",  # Waldemath Moon - hypothetical second moon of Earth
    "calc_planet_x_pickering",  # Pickering's Planet O/X prediction (1919)
    "calc_white_moon_position",  # White Moon (Selena) - opposite to Black Moon Lilith
    "calc_uranian_planet",
    "calc_hypothetical_position",
    "URANIAN_KEPLERIAN_ELEMENTS",
    "TRANSPLUTO_KEPLERIAN_ELEMENTS",
    "VULCAN_ELEMENTS",
    "WALDEMATH_ELEMENTS",  # Waldemath Moon orbital elements
    "PICKERING_PLANET_X_ELEMENTS",  # Pickering's Planet O/X orbital elements
    # seorbel.txt parser for custom hypothetical bodies
    "parse_seorbel",
    "get_bundled_seorbel_path",
    "load_bundled_seorbel",
    "SeorbelElements",
    "TPolynomial",
    "get_seorbel_body_by_name",
    "calc_seorbel_position",
    # Planetary moons (Galilean moons, Titan, etc.)
    "register_moon_spk",
    "unregister_moon_spk",
    "list_registered_moons",
    "get_moon_name",
    "is_planetary_moon",
    "get_moon_coverage",
    "calc_moon_position",
    "close_moon_kernels",
    "MOON_NAIF_MAP",
    "MOON_NAMES",
    "MOON_PARENT_MAP",
    # REBOUND/ASSIST n-body integration
    "ReboundIntegrator",
    "AssistEphemConfig",
    "PropagationResult",
    "check_rebound_available",
    "check_assist_available",
    "get_rebound_version",
    "get_assist_version",
    "propagate_orbit_rebound",
    "propagate_orbit_assist",
    "propagate_trajectory",
    "compare_with_keplerian",
]
