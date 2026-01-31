"""
REBOUND/ASSIST integration for ephemeris-quality asteroid orbit propagation.

This module provides high-precision n-body gravitational integration for asteroid
and comet orbit propagation using the REBOUND and ASSIST libraries.

REBOUND is a multi-purpose N-body integrator with several integrator choices:
- IAS15: High-accuracy adaptive integrator (15th order Gauss-Radau)
- WHFast: Fast symplectic Wisdom-Holman integrator
- MERCURIUS: Hybrid symplectic integrator for close encounters
- TRACE: Hybrid reversible integrator for arbitrary close encounters

ASSIST extends REBOUND for ephemeris-quality integrations of Solar System test
particles, matching JPL's small body integrator precision:
- Uses JPL DE440/DE441 ephemeris for planet positions
- Includes Sun, Moon, 8 planets, and 16 massive asteroids
- Gravitational harmonics (J2, J3, J4) for Earth and Sun
- General relativistic corrections
- Non-gravitational effects (Marsden model)
- Verified to machine precision against JPL Horizons

INSTALLATION:
    pip install rebound          # Base n-body integrator (~3 MB)
    pip install assist           # Ephemeris-quality extension (~3 MB)

    # ASSIST requires ephemeris data files (~1 GB total):
    mkdir -p data
    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440
    curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o data/sb441-n16.bsp

USAGE:
    >>> from libephemeris.rebound_integration import (
    ...     propagate_orbit_assist,
    ...     propagate_orbit_rebound,
    ...     ReboundIntegrator,
    ... )
    >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, SE_CERES
    >>>
    >>> # Propagate Ceres orbit using ASSIST (ephemeris-quality)
    >>> elements = MINOR_BODY_ELEMENTS[SE_CERES]
    >>> jd_start = 2460000.5  # Starting JD
    >>> jd_end = 2460365.5    # One year later
    >>> x, y, z, vx, vy, vz = propagate_orbit_assist(elements, jd_start, jd_end)

PRECISION:
    - ASSIST: Sub-arcsecond accuracy, matching JPL Horizons
    - REBOUND IAS15: ~0.001-1 arcsecond (depends on timestep and duration)
    - REBOUND WHFast: ~1-10 arcseconds (better energy conservation for long term)
    - Keplerian (minor_bodies.py): ~10-30 arcseconds for main belt asteroids

REFERENCES:
    - REBOUND: https://rebound.readthedocs.io/
    - ASSIST: https://assist.readthedocs.io/
    - Rein & Liu 2012, A&A 537, A128 (REBOUND)
    - Rein & Spiegel 2015, MNRAS 446, 1424 (IAS15)
    - Rein & Tamayo 2015, MNRAS 452, 376 (WHFast)
    - Holman et al. 2023 (ASSIST)
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Tuple, List, Callable

if TYPE_CHECKING:
    from .minor_bodies import OrbitalElements


class ReboundIntegrator(Enum):
    """Available integrators in REBOUND.

    Each integrator has different trade-offs between speed, accuracy, and
    use cases.
    """

    IAS15 = "ias15"
    """15th order Gauss-Radau integrator with adaptive timestep.
    
    Best for: High-accuracy integrations, close encounters, eccentric orbits.
    Accuracy: Machine precision for most problems.
    Speed: Moderate (adaptive, slows for close encounters).
    """

    WHFAST = "whfast"
    """Fast Wisdom-Holman symplectic integrator.
    
    Best for: Long-term integrations of planetary systems.
    Accuracy: O(dt^2) per step, but bounded errors (symplectic).
    Speed: Fast (fixed timestep).
    Note: Not suitable for close encounters.
    """

    MERCURIUS = "mercurius"
    """Hybrid symplectic integrator for close encounters.
    
    Best for: Systems with occasional close encounters.
    Accuracy: Symplectic most of time, IAS15 during encounters.
    Speed: Fast until close encounters occur.
    """

    TRACE = "trace"
    """Hybrid TRAnsits and Close Encounters integrator.
    
    Best for: Arbitrary close encounters, reversible integration.
    Accuracy: High accuracy even with many close encounters.
    Speed: Moderate.
    """

    LEAPFROG = "leapfrog"
    """Standard leapfrog/Verlet integrator.
    
    Best for: Simple problems, educational use.
    Accuracy: O(dt^2) per step.
    Speed: Fast.
    """


@dataclass
class AssistEphemConfig:
    """Configuration for ASSIST ephemeris data files.

    ASSIST requires two ephemeris data files:
    1. Planet ephemeris (DE440/DE441): ~400 MB - 2.6 GB
    2. Asteroid perturber ephemeris: ~600 MB

    These files can be downloaded from:
    - https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/
    - https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/

    Attributes:
        planets_file: Path to planet ephemeris file (e.g., linux_p1550p2650.440)
        asteroids_file: Path to asteroid perturber file (e.g., sb441-n16.bsp)
        data_dir: Directory containing ephemeris files (alternative to full paths)
    """

    planets_file: Optional[str] = None
    asteroids_file: Optional[str] = None
    data_dir: Optional[str] = None

    def __post_init__(self):
        """Resolve file paths from environment or defaults."""
        # Try environment variable
        env_dir = os.environ.get("ASSIST_DIR")

        if self.data_dir is None and env_dir:
            self.data_dir = env_dir

        # Set default filenames if not provided
        if self.data_dir:
            base = Path(self.data_dir)
            if self.planets_file is None:
                # Try different planet ephemeris files in order of preference
                for fname in [
                    "de441.bsp",
                    "de440.bsp",
                    "linux_m13000p17000.441",
                    "linux_p1550p2650.440",
                ]:
                    if (base / fname).exists():
                        self.planets_file = str(base / fname)
                        break

            if self.asteroids_file is None:
                for fname in ["sb441-n16.bsp"]:
                    if (base / fname).exists():
                        self.asteroids_file = str(base / fname)
                        break


@dataclass
class PropagationResult:
    """Result of orbit propagation.

    Contains position and velocity vectors in heliocentric ecliptic J2000
    coordinates at the target time.

    Attributes:
        x, y, z: Position components in AU
        vx, vy, vz: Velocity components in AU/day
        jd_tt: Julian Date (TT) of the result
        ecliptic_lon: Ecliptic longitude in degrees
        ecliptic_lat: Ecliptic latitude in degrees
        distance: Heliocentric distance in AU
    """

    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    jd_tt: float

    @property
    def ecliptic_lon(self) -> float:
        """Ecliptic longitude in degrees (0-360)."""
        return math.degrees(math.atan2(self.y, self.x)) % 360.0

    @property
    def ecliptic_lat(self) -> float:
        """Ecliptic latitude in degrees (-90 to +90)."""
        return math.degrees(math.asin(self.z / self.distance))

    @property
    def distance(self) -> float:
        """Heliocentric distance in AU."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def to_tuple(self) -> Tuple[float, float, float, float, float, float]:
        """Return (x, y, z, vx, vy, vz) tuple."""
        return (self.x, self.y, self.z, self.vx, self.vy, self.vz)


def check_rebound_available() -> bool:
    """Check if REBOUND is installed and available.

    Returns:
        bool: True if REBOUND can be imported, False otherwise.
    """
    try:
        import rebound

        return True
    except ImportError:
        return False


def check_assist_available() -> bool:
    """Check if ASSIST is installed and available.

    Returns:
        bool: True if ASSIST can be imported, False otherwise.
    """
    try:
        import assist

        return True
    except ImportError:
        return False


def get_rebound_version() -> Optional[str]:
    """Get the installed REBOUND version.

    Returns:
        str: Version string if installed, None otherwise.
    """
    try:
        import rebound

        return rebound.__version__
    except ImportError:
        return None


def get_assist_version() -> Optional[str]:
    """Get the installed ASSIST version.

    Returns:
        str: Version string if installed, None otherwise.
    """
    try:
        import assist

        return assist.__version__
    except ImportError:
        return None


def elements_to_rebound_orbit(
    elements: "OrbitalElements",
    jd_tt: float,
) -> dict:
    """Convert OrbitalElements to REBOUND orbital parameters.

    REBOUND uses different conventions for orbital elements:
    - Semi-major axis in simulation units (typically AU)
    - Angles in radians (we convert from degrees)
    - Mean anomaly or true anomaly

    Args:
        elements: OrbitalElements from minor_bodies module
        jd_tt: Julian Date (TT) at which to evaluate elements

    Returns:
        dict: Parameters suitable for rebound.Particle.add(...)
    """
    from .minor_bodies import apply_secular_perturbations

    # Get perturbed elements at target time
    omega, Omega, M, n = apply_secular_perturbations(
        elements, jd_tt, include_perturbations=True
    )

    return {
        "a": elements.a,  # Semi-major axis in AU
        "e": elements.e,  # Eccentricity
        "inc": math.radians(elements.i),  # Inclination in radians
        "omega": math.radians(omega),  # Argument of perihelion in radians
        "Omega": math.radians(Omega),  # Longitude of ascending node in radians
        "M": math.radians(M),  # Mean anomaly in radians
    }


def propagate_orbit_rebound(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    integrator: ReboundIntegrator = ReboundIntegrator.IAS15,
    dt: Optional[float] = None,
) -> PropagationResult:
    """Propagate an orbit using REBOUND N-body integration.

    Uses REBOUND to integrate the orbit in the gravitational field of the Sun
    only (2-body problem). For ephemeris-quality precision including planetary
    perturbations, use propagate_orbit_assist() instead.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        integrator: Which REBOUND integrator to use (default: IAS15)
        dt: Timestep in days for fixed-step integrators (default: auto)

    Returns:
        PropagationResult: Position and velocity at jd_end

    Raises:
        ImportError: If REBOUND is not installed
        ValueError: If elements are invalid

    Example:
        >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, SE_CERES
        >>> elements = MINOR_BODY_ELEMENTS[SE_CERES]
        >>> result = propagate_orbit_rebound(elements, 2460000.5, 2460365.5)
        >>> print(f"Ceres at {result.ecliptic_lon:.4f} deg, {result.distance:.4f} AU")
    """
    try:
        import rebound
    except ImportError:
        raise ImportError(
            "REBOUND is required for n-body orbit propagation. "
            "Install with: pip install rebound"
        )

    # Gaussian gravitational constant squared (k^2)
    # This gives G*M_sun in AU^3/day^2
    k_squared = 0.00029591220828559  # AU^3/day^2

    # Create simulation
    sim = rebound.Simulation()
    sim.units = ("day", "AU", "Msun")
    sim.G = k_squared  # GM_sun = k^2 for Msun = 1

    # Set integrator
    sim.integrator = integrator.value

    # Set timestep for fixed-step integrators
    if dt is not None:
        sim.dt = dt
    elif integrator in (ReboundIntegrator.WHFAST, ReboundIntegrator.LEAPFROG):
        # Default to ~1/20 of orbital period
        period = 365.25 * elements.a**1.5  # Approximate period in days
        sim.dt = period / 20.0

    # Add the Sun at origin
    sim.add(m=1.0, hash="Sun")

    # Convert orbital elements to REBOUND format
    orb_params = elements_to_rebound_orbit(elements, jd_start)

    # Add test particle with orbital elements
    sim.add(
        m=0.0,  # Test particle (no mass)
        primary=sim.particles["Sun"],
        hash=elements.name,
        **orb_params,
    )

    # Move to heliocentric frame (Sun at origin)
    sim.move_to_hel()

    # Set initial time (relative, simulation starts at t=0)
    sim.t = 0.0

    # Integration time in days
    integration_time = jd_end - jd_start

    # Integrate
    sim.integrate(integration_time)

    # Get particle state
    p = sim.particles[1]  # The asteroid (index 1, after Sun)

    return PropagationResult(
        x=p.x,
        y=p.y,
        z=p.z,
        vx=p.vx,
        vy=p.vy,
        vz=p.vz,
        jd_tt=jd_end,
    )


def propagate_orbit_assist(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    ephem_config: Optional[AssistEphemConfig] = None,
    include_non_gravitational: bool = False,
    A1: float = 0.0,
    A2: float = 0.0,
    A3: float = 0.0,
) -> PropagationResult:
    """Propagate an orbit using ASSIST for ephemeris-quality integration.

    ASSIST provides the highest precision orbit propagation by including:
    - Sun, Moon, 8 planets from JPL DE440/DE441 ephemeris
    - 16 massive asteroids (Ceres, Vesta, Pallas, etc.)
    - Earth and Sun gravitational harmonics (J2, J3, J4)
    - General relativistic corrections
    - Optional non-gravitational forces (Marsden model)

    The precision matches JPL's small body integrator to machine precision.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        ephem_config: Configuration for ephemeris data files
        include_non_gravitational: Include Marsden non-gravitational forces
        A1, A2, A3: Marsden parameters for non-gravitational acceleration
            (only used if include_non_gravitational=True)

    Returns:
        PropagationResult: Position and velocity at jd_end

    Raises:
        ImportError: If ASSIST or REBOUND is not installed
        FileNotFoundError: If ephemeris files are not found
        ValueError: If elements are invalid

    Example:
        >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, SE_CERES
        >>> from libephemeris.rebound_integration import propagate_orbit_assist
        >>> elements = MINOR_BODY_ELEMENTS[SE_CERES]
        >>> result = propagate_orbit_assist(elements, 2460000.5, 2460365.5)
        >>> print(f"Ceres: lon={result.ecliptic_lon:.6f} deg")

    Note:
        Requires ~1 GB of ephemeris data files. See module docstring for
        download instructions.
    """
    try:
        import rebound
    except ImportError:
        raise ImportError(
            "REBOUND is required for ASSIST. Install with: pip install rebound"
        )

    try:
        import assist
    except ImportError:
        raise ImportError(
            "ASSIST is required for ephemeris-quality integration. "
            "Install with: pip install assist\n"
            "See also: https://assist.readthedocs.io/"
        )

    # Configure ephemeris paths
    if ephem_config is None:
        ephem_config = AssistEphemConfig()

    # Verify ephemeris files exist
    if ephem_config.planets_file is None:
        raise FileNotFoundError(
            "Planet ephemeris file not found. Download from:\n"
            "curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440\n"
            "Set ASSIST_DIR environment variable or pass ephem_config."
        )

    if not os.path.exists(ephem_config.planets_file):
        raise FileNotFoundError(
            f"Planet ephemeris file not found: {ephem_config.planets_file}"
        )

    if ephem_config.asteroids_file and not os.path.exists(ephem_config.asteroids_file):
        raise FileNotFoundError(
            f"Asteroid ephemeris file not found: {ephem_config.asteroids_file}"
        )

    # Create ASSIST ephemeris object
    if ephem_config.asteroids_file:
        ephem = assist.Ephem(ephem_config.planets_file, ephem_config.asteroids_file)
    else:
        ephem = assist.Ephem(ephem_config.planets_file)

    # Create REBOUND simulation
    sim = rebound.Simulation()

    # Attach ASSIST extras
    extras = assist.Extras(sim, ephem)

    # Configure non-gravitational forces if requested
    if include_non_gravitational:
        # Set Marsden model parameters
        # A1: radial, A2: transverse, A3: normal components
        extras.particle_params = {
            "A1": A1,
            "A2": A2,
            "A3": A3,
        }

    # Convert orbital elements to Cartesian state
    orb_params = elements_to_rebound_orbit(elements, jd_start)

    # Add test particle
    sim.add(
        m=0.0,  # Test particle
        **orb_params,
    )

    # Set initial time (JD - reference JD)
    sim.t = jd_start - ephem.jd_ref

    # Integration time
    t_end = jd_end - ephem.jd_ref

    # Integrate
    sim.integrate(t_end)

    # Get particle state
    p = sim.particles[0]

    return PropagationResult(
        x=p.x,
        y=p.y,
        z=p.z,
        vx=p.vx,
        vy=p.vy,
        vz=p.vz,
        jd_tt=jd_end,
    )


def propagate_trajectory(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    num_points: int = 100,
    use_assist: bool = True,
    ephem_config: Optional[AssistEphemConfig] = None,
    integrator: ReboundIntegrator = ReboundIntegrator.IAS15,
) -> List[PropagationResult]:
    """Propagate an orbit and return positions at multiple times.

    Useful for generating ephemerides or plotting orbital trajectories.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        num_points: Number of output points (evenly spaced in time)
        use_assist: Use ASSIST for ephemeris-quality (if available)
        ephem_config: ASSIST ephemeris configuration
        integrator: REBOUND integrator (only used if use_assist=False)

    Returns:
        List[PropagationResult]: Positions and velocities at each time

    Example:
        >>> trajectory = propagate_trajectory(elements, jd_start, jd_end, num_points=365)
        >>> for point in trajectory:
        ...     print(f"JD {point.jd_tt}: lon={point.ecliptic_lon:.4f}")
    """
    results = []

    # Time step between output points
    dt = (jd_end - jd_start) / (num_points - 1) if num_points > 1 else 0

    if use_assist and check_assist_available():
        try:
            import rebound
            import assist

            # Configure ephemeris
            if ephem_config is None:
                ephem_config = AssistEphemConfig()

            if ephem_config.planets_file and os.path.exists(ephem_config.planets_file):
                # Use ASSIST
                if ephem_config.asteroids_file:
                    ephem = assist.Ephem(
                        ephem_config.planets_file, ephem_config.asteroids_file
                    )
                else:
                    ephem = assist.Ephem(ephem_config.planets_file)

                sim = rebound.Simulation()
                extras = assist.Extras(sim, ephem)

                orb_params = elements_to_rebound_orbit(elements, jd_start)
                sim.add(m=0.0, **orb_params)

                sim.t = jd_start - ephem.jd_ref

                for i in range(num_points):
                    t = jd_start + i * dt
                    sim.integrate(t - ephem.jd_ref)
                    p = sim.particles[0]
                    results.append(
                        PropagationResult(
                            x=p.x, y=p.y, z=p.z, vx=p.vx, vy=p.vy, vz=p.vz, jd_tt=t
                        )
                    )

                return results
        except (ImportError, FileNotFoundError):
            pass  # Fall back to REBOUND

    # Use REBOUND without ASSIST
    if not check_rebound_available():
        raise ImportError(
            "REBOUND is required for orbit propagation. "
            "Install with: pip install rebound"
        )

    import rebound

    k_squared = 0.00029591220828559
    sim = rebound.Simulation()
    sim.units = ("day", "AU", "Msun")
    sim.G = k_squared
    sim.integrator = integrator.value

    sim.add(m=1.0, hash="Sun")

    orb_params = elements_to_rebound_orbit(elements, jd_start)
    sim.add(m=0.0, primary=sim.particles["Sun"], **orb_params)
    sim.move_to_hel()
    sim.t = 0.0

    for i in range(num_points):
        t = i * dt
        sim.integrate(t)
        p = sim.particles[1]
        results.append(
            PropagationResult(
                x=p.x, y=p.y, z=p.z, vx=p.vx, vy=p.vy, vz=p.vz, jd_tt=jd_start + t
            )
        )

    return results


def compare_with_keplerian(
    elements: "OrbitalElements",
    jd_tt: float,
    use_assist: bool = False,
    ephem_config: Optional[AssistEphemConfig] = None,
) -> dict:
    """Compare REBOUND/ASSIST integration with Keplerian approximation.

    Useful for evaluating the improvement in precision from using numerical
    integration vs. the analytical Keplerian method with secular perturbations.

    Args:
        elements: OrbitalElements to propagate
        jd_tt: Target Julian Date (TT)
        use_assist: Use ASSIST if available (else REBOUND only)
        ephem_config: ASSIST ephemeris configuration

    Returns:
        dict: Comparison results including positions and differences

    Example:
        >>> comparison = compare_with_keplerian(elements, jd_tt)
        >>> print(f"Position difference: {comparison['angular_sep_arcsec']:.2f} arcsec")
    """
    from .minor_bodies import calc_minor_body_position

    # Keplerian result
    x_kep, y_kep, z_kep = calc_minor_body_position(elements, jd_tt)
    r_kep = math.sqrt(x_kep**2 + y_kep**2 + z_kep**2)
    lon_kep = math.degrees(math.atan2(y_kep, x_kep)) % 360.0
    lat_kep = math.degrees(math.asin(z_kep / r_kep))

    # REBOUND/ASSIST result
    jd_start = elements.epoch

    if use_assist and check_assist_available():
        try:
            result = propagate_orbit_assist(elements, jd_start, jd_tt, ephem_config)
        except (ImportError, FileNotFoundError):
            result = propagate_orbit_rebound(elements, jd_start, jd_tt)
    else:
        result = propagate_orbit_rebound(elements, jd_start, jd_tt)

    # Angular separation
    d_lon = result.ecliptic_lon - lon_kep
    if d_lon > 180:
        d_lon -= 360
    elif d_lon < -180:
        d_lon += 360

    d_lat = result.ecliptic_lat - lat_kep

    # Approximate angular separation in arcseconds
    cos_lat = math.cos(math.radians((lat_kep + result.ecliptic_lat) / 2))
    angular_sep = math.sqrt((d_lon * cos_lat) ** 2 + d_lat**2) * 3600  # arcsec

    return {
        "keplerian": {
            "x": x_kep,
            "y": y_kep,
            "z": z_kep,
            "lon": lon_kep,
            "lat": lat_kep,
            "dist": r_kep,
        },
        "nbody": {
            "x": result.x,
            "y": result.y,
            "z": result.z,
            "lon": result.ecliptic_lon,
            "lat": result.ecliptic_lat,
            "dist": result.distance,
        },
        "delta_lon_arcsec": d_lon * 3600,
        "delta_lat_arcsec": d_lat * 3600,
        "angular_sep_arcsec": angular_sep,
        "delta_dist_au": result.distance - r_kep,
        "method": "assist" if use_assist else "rebound",
        "propagation_days": jd_tt - jd_start,
    }
