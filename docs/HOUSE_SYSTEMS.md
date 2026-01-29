# House System Algorithms and Mathematical Formulas

This document provides comprehensive mathematical documentation for all 19 house systems implemented in LibEphemeris. Each section includes the underlying algorithm, mathematical formulas, and implementation notes.

## Table of Contents

1. [Terminology and Definitions](#terminology-and-definitions)
2. [Placidus (P)](#placidus-p)
3. [Koch (K)](#koch-k)
4. [Regiomontanus (R)](#regiomontanus-r)
5. [Campanus (C)](#campanus-c)
6. [Porphyry (O)](#porphyry-o)
7. [Equal (A/E)](#equal-ae)
8. [Whole Sign (W)](#whole-sign-w)
9. [Alcabitius (B)](#alcabitius-b)
10. [Topocentric/Polich-Page (T)](#topocentricpolich-page-t)
11. [Morinus (M)](#morinus-m)
12. [Meridian/Axial Rotation (X)](#meridianaxial-rotation-x)
13. [Vehlow (V)](#vehlow-v)
14. [Carter Poli-Equatorial (F)](#carter-poli-equatorial-f)
15. [Gauquelin Sectors (G)](#gauquelin-sectors-g)
16. [Horizontal/Azimuthal (H)](#horizontalazimuthal-h)
17. [Krusinski-Pisa (U)](#krusinski-pisa-u)
18. [Natural Gradient (N)](#natural-gradient-n)
19. [APC (Y)](#apc-y)

---

## Terminology and Definitions

### Coordinate Systems

- **Ecliptic coordinates (λ, β)**: Longitude and latitude measured along the ecliptic plane
- **Equatorial coordinates (α, δ)**: Right Ascension (RA) and Declination
- **Horizontal coordinates (A, h)**: Azimuth and Altitude

### Key Variables

| Symbol | Definition | Units |
|--------|-----------|-------|
| λ | Ecliptic longitude | degrees (0-360) |
| ε | Obliquity of the ecliptic | degrees (~23.44°) |
| φ | Geographic latitude | degrees (-90 to +90) |
| θ | Local Sidereal Time (LST) | degrees |
| ARMC | Right Ascension of MC | degrees |
| α | Right Ascension | degrees |
| δ | Declination | degrees |
| AD | Ascensional Difference | degrees |
| SA | Semi-diurnal Arc | degrees |
| OA | Oblique Ascension | degrees |
| H | Hour Angle | degrees |

### Fundamental Formulas

**Ecliptic to Equatorial Conversion:**
```
tan(α) = sin(λ) · cos(ε) / cos(λ)
sin(δ) = sin(λ) · sin(ε)
```

**Equatorial to Ecliptic Conversion:**
```
tan(λ) = sin(α) / (cos(α) · cos(ε))
```

**Ascensional Difference:**
```
sin(AD) = tan(φ) · tan(δ)
```

**Semi-diurnal Arc:**
```
SA = 90° + AD  (diurnal, above horizon)
SA = 90° - AD  (nocturnal, below horizon)
```

---

## Placidus (P)

**Historical Origin:** Placidus de Titis (1603-1668), Italian mathematician and monk.

**Conceptual Basis:** Time-based division of the semi-diurnal and semi-nocturnal arcs into thirds.

### Algorithm

1. Calculate Semi-diurnal Arc (SA) for each point
2. Divide SA into three equal time segments
3. Iteratively find ecliptic longitude at each time division
4. Houses 11, 12 from MC to Asc (above horizon)
5. Houses 2, 3 from Asc to IC (below horizon)
6. Opposite houses (5, 6, 8, 9) are 180° from (11, 12, 2, 3)

### Mathematical Formulas

**House 11 (1/3 of semi-arc from MC):**
```
H₁₁ = SA/3 = (90° + AD)/3
RA₁₁ = ARMC + H₁₁
```

**House 12 (2/3 of semi-arc from MC):**
```
H₁₂ = 2·SA/3 = 2·(90° + AD)/3
RA₁₂ = ARMC + H₁₂
```

**House 2 (2/3 of nocturnal semi-arc from IC):**
```
H₂ = 2·(90° - AD)/3
RA₂ = ARMC + 180° - H₂
```

**House 3 (1/3 of nocturnal semi-arc from IC):**
```
H₃ = (90° - AD)/3
RA₃ = ARMC + 180° - H₃
```

**Iterative Solution:**
For each cusp, iterate until convergence (threshold: 1e-7°):
```
1. tan(δ) = sin(RA) · tan(ε)
2. AD = arcsin(tan(φ) · tan(δ))
3. RA_new = ARMC + H (computed from AD)
4. If |RA_new - RA| < 1e-7°, stop
```

**RA to Ecliptic Longitude:**
```
tan(λ) = sin(RA) / (cos(RA) · cos(ε))
λ = atan2(sin(RA), cos(RA) · cos(ε))
```

### Polar Limitation

Fails when |φ| + ε > 90° because tan(φ) · tan(δ) > 1, making AD undefined.

---

## Koch (K)

**Historical Origin:** Walter Koch (1895-1970), German astrologer. Also called GOH (Geburtsort-Häusersystem) or Birthplace system.

**Conceptual Basis:** Trisection of Oblique Ascension intervals between angles.

### Algorithm

1. Calculate Oblique Ascension for MC, Asc, and IC
2. Divide OA intervals into thirds
3. Iteratively solve for ecliptic longitude at each target OA
4. Calculate opposite houses by adding 180°

### Mathematical Formulas

**Oblique Ascension:**
```
OA = RA - AD
where AD = arcsin(tan(φ) · tan(δ))
```

**For houses 11, 12 (MC to Asc quadrant):**
```
OA₁₁ = OA_MC + (OA_Asc - OA_MC)/3
OA₁₂ = OA_MC + 2·(OA_Asc - OA_MC)/3
```

**For houses 2, 3 (Asc to IC quadrant):**
```
OA₂ = OA_Asc + (OA_IC - OA_Asc)/3
OA₃ = OA_Asc + 2·(OA_IC - OA_Asc)/3
```

**Iterative Solution for target OA:**
```
1. RA = target_OA (initial guess)
2. tan(δ) = sin(RA) · tan(ε)
3. AD = arcsin(tan(φ) · tan(δ))
4. RA_new = target_OA + AD
5. Repeat until |RA_new - RA| < 1e-7°
6. Convert final RA to ecliptic longitude
```

### Polar Limitation

Same as Placidus: fails when |φ| + ε > 90°.

---

## Regiomontanus (R)

**Historical Origin:** Johannes Müller von Königsberg (Regiomontanus, 1436-1476).

**Conceptual Basis:** Equal 30° divisions of the celestial equator, projected to ecliptic via great circles through the north and south points of the horizon.

### Algorithm

1. Divide equator into 30° segments from ARMC
2. Calculate pole of projection for each segment
3. Project equator point to ecliptic using spherical trigonometry

### Mathematical Formulas

**Pole of projection:**
```
tan(P) = tan(φ) · sin(H)
where H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3
```

**Right Ascension offset:**
```
R = ARMC + H - 90°
```

**Ecliptic longitude:**
```
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Campanus (C)

**Historical Origin:** Giovanni di Campani (13th century).

**Conceptual Basis:** Equal 30° divisions of the prime vertical, projected onto the ecliptic.

### Algorithm

1. Divide prime vertical into 30° segments
2. Transform each segment to equatorial hour angle
3. Apply Regiomontanus-style projection

### Mathematical Formulas

**Prime vertical to hour angle transformation:**
```
tan(H_eff) = tan(h_pv) · cos(φ)
where h_pv = 30°, 60°, 120°, 150° (prime vertical offset)
```

**Pole calculation:**
```
tan(P) = tan(φ) · sin(H_eff)
```

**Ecliptic projection:**
```
R = ARMC + H_eff - 90°
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Porphyry (O)

**Historical Origin:** Porphyry of Tyre (234-305 CE), Neoplatonic philosopher.

**Conceptual Basis:** Space-based trisection of ecliptic arc between angles.

### Algorithm

1. Divide ecliptic arc from MC to Asc into three equal parts
2. Divide ecliptic arc from Asc to IC into three equal parts
3. Opposite houses are 180° apart

### Mathematical Formulas

**Quadrant 10→1 (MC to Asc):**
```
step = (λ_Asc - λ_MC) mod 360° / 3
λ₁₁ = λ_MC + step
λ₁₂ = λ_MC + 2·step
```

**Quadrant 1→4 (Asc to IC):**
```
step = (λ_IC - λ_Asc) mod 360° / 3
λ₂ = λ_Asc + step
λ₃ = λ_Asc + 2·step
```

**Opposite houses:**
```
λ₅ = (λ₁₁ + 180°) mod 360°
λ₆ = (λ₁₂ + 180°) mod 360°
λ₈ = (λ₂ + 180°) mod 360°
λ₉ = (λ₃ + 180°) mod 360°
```

**Advantages:** Works at all latitudes, computationally simple.

---

## Equal (A/E)

**Historical Origin:** Ancient Hellenistic astrology.

**Conceptual Basis:** Each house is exactly 30° of ecliptic longitude.

### Mathematical Formula

```
λᵢ = (λ_Asc + (i-1)·30°) mod 360°
where i = 1, 2, ..., 12
```

**Properties:**
- Simplest calculation
- Works at all latitudes
- MC may not coincide with 10th house cusp

---

## Whole Sign (W)

**Historical Origin:** Oldest house system, used in Hellenistic, Indian, and Medieval astrology.

**Conceptual Basis:** Each house is one complete zodiac sign; House 1 is the sign containing the Ascendant.

### Mathematical Formula

```
start = floor(λ_Asc / 30°) × 30°
λᵢ = (start + (i-1)·30°) mod 360°
where i = 1, 2, ..., 12
```

**Example:** If Asc = 15° Taurus (45°), then:
- House 1 starts at 30° (0° Taurus)
- House 2 starts at 60° (0° Gemini)
- etc.

---

## Alcabitius (B)

**Historical Origin:** Abd al-Aziz al-Qabisi (Alcabitius, 10th century), Arab astrologer.

**Conceptual Basis:** Time trisection of Ascendant's diurnal arc, projected by hour circles.

### Algorithm

1. Calculate Right Ascension of Ascendant
2. Divide RA intervals between angles
3. Convert each RA division back to ecliptic longitude

### Mathematical Formulas

**Quadrant MC to Asc:**
```
arc = (RA_Asc - ARMC) mod 360°
step = arc / 3
RA₁₁ = ARMC + step
RA₁₂ = ARMC + 2·step
```

**Quadrant Asc to IC:**
```
arc = (RA_IC - RA_Asc) mod 360°
step = arc / 3
RA₂ = RA_Asc + step
RA₃ = RA_Asc + 2·step
```

**RA to ecliptic conversion:**
```
λ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Topocentric/Polich-Page (T)

**Historical Origin:** Wendel Polich and A.P. Nelson Page (1961).

**Conceptual Basis:** Modified pole method to account for observer's actual position on Earth's surface.

### Algorithm

Similar to Regiomontanus but with modified pole factors.

### Mathematical Formulas

**Pole factors:**
```
For houses 11, 3: factor = 1/3
For houses 12, 2: factor = 2/3
```

**Pole calculation:**
```
tan(P) = tan(φ) · factor
```

**Ecliptic projection:**
```
R = ARMC + H - 90°
where H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Morinus (M)

**Historical Origin:** Jean-Baptiste Morin (1583-1656), French astrologer.

**Conceptual Basis:** Equal 30° divisions on celestial equator, projected to ecliptic via ecliptic poles.

### Algorithm

Project equator points (ARMC + n·30°) to ecliptic.

### Mathematical Formula

```
For each house cusp at RA = ARMC + (i-1)·30°:
tan(λ) = tan(RA) · cos(ε)
λ = atan2(sin(RA)·cos(ε), cos(RA))
```

**Properties:**
- Location-independent (ignores latitude)
- Morinus 10th cusp differs from standard MC

---

## Meridian/Axial Rotation (X)

**Historical Origin:** Also known as Zariel system.

**Conceptual Basis:** Equal 30° RA divisions from MC, projected to ecliptic via celestial poles.

### Mathematical Formula

```
For each house cusp at RA = ARMC + (i-10)·30°:
λ = atan2(sin(RA), cos(RA)·cos(ε))
```

**Difference from Morinus:**
```
Morinus: λ = atan2(sin(RA)·cos(ε), cos(RA))
Meridian: λ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Vehlow (V)

**Historical Origin:** Johannes Vehlow (1890-1958), German astrologer.

**Conceptual Basis:** Equal houses with Ascendant at center (15°) of House 1 rather than at cusp.

### Mathematical Formula

```
start = (λ_Asc - 15°) mod 360°
λᵢ = (start + (i-1)·30°) mod 360°
```

**Result:** Ascendant falls at 15° into House 1.

---

## Carter Poli-Equatorial (F)

**Historical Origin:** Charles Carter (1887-1968), English astrologer.

**Conceptual Basis:** Equal 30° divisions on equator starting from RA of Ascendant, projected to ecliptic.

### Mathematical Formula

```
RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

For each house i (1-12):
RA = RA_Asc + (i-1)·30°
λᵢ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Gauquelin Sectors (G)

**Historical Origin:** Michel and Françoise Gauquelin (20th century), for statistical research.

**Conceptual Basis:** 36 sectors based on true diurnal and nocturnal arc division.

### Algorithm

1. Calculate diurnal arc (Asc → MC → Desc): 18 sectors
2. Calculate nocturnal arc (Desc → IC → Asc): 18 sectors
3. Map 36 sectors to 12 houses (3 sectors per house)
4. Use middle sector as house cusp

### Mathematical Formulas

**Arc calculation:**
```
arc_diurnal = (RA_Desc - RA_Asc) mod 360°
arc_nocturnal = 360° - arc_diurnal
```

**Sector RA positions (1-36):**
```
For sectors 1-18 (diurnal):
  RA_sector = RA_Asc + (sector-1)/18 × arc_diurnal

For sectors 19-36 (nocturnal):
  RA_sector = RA_Desc + (sector-19)/18 × arc_nocturnal
```

**House cusp mapping:**
```
House 1 cusp = sector 2 longitude
House 2 cusp = sector 5 longitude
House n cusp = sector (3n-1) longitude
```

**RA to ecliptic:**
```
λ = atan2(sin(RA), cos(RA)·cos(ε))
```

### Polar Limitation

Fails within polar circle due to undefined diurnal/nocturnal arcs.

---

## Horizontal/Azimuthal (H)

**Historical Origin:** Traditional horizon-based system.

**Conceptual Basis:** House circles based on horizon plane divisions.

### Algorithm

Uses co-latitude transformation and Campanus-like projection:

1. Transform latitude to co-latitude: `φ' = 90° - φ`
2. Calculate intermediate azimuths
3. Apply ascendant formula with modified angles

### Mathematical Formulas

**Co-latitude:**
```
φ' = 90° - φ  (for φ > 0)
φ' = -90° - φ (for φ < 0)
```

**Intermediate azimuths:**
```
fh1 = arcsin(sin(φ')/2)
fh2 = arcsin(√3/2 · sin(φ'))
xh1 = arctan(√3 / cos(φ'))
xh2 = arctan(1/(√3 · cos(φ')))
```

**House cusps (using modified ARMC):**
```
θ' = ARMC + 180°
cusps[11] = Asc(θ' + 90° - xh1, ε, φ, fh1)
cusps[12] = Asc(θ' + 90° - xh2, ε, φ, fh2)
cusps[1]  = Asc(θ' + 90°, ε, φ, φ')
cusps[2]  = Asc(θ' + 90° + xh2, ε, φ, fh2)
cusps[3]  = Asc(θ' + 90° + xh1, ε, φ, fh1)
```

---

## Krusinski-Pisa (U)

**Historical Origin:** Bogdan Krusinski (modern Polish astrologer).

**Conceptual Basis:** Great circle passing through Ascendant and Zenith, divided into 12 equal parts.

### Algorithm

1. Transform Ascendant from ecliptic to equatorial coordinates
2. Rotate to align with meridian
3. Transform to horizontal coordinates
4. Create Asc-Zenith great circle
5. Divide into 30° segments
6. Transform each segment back to ecliptic

### Mathematical Formulas

**Coordinate rotation (cotrans):**
```
x' = x
y' = y·cos(θ) + z·sin(θ)
z' = -y·sin(θ) + z·cos(θ)
```

**Algorithm steps:**
```
A0. x = [λ_Asc, 0, 1]  (ecliptic coords)
A1. x = cotrans(x, -ε)  (to equatorial)
A2. x[0] = x[0] - (ARMC - 90°)  (rotate)
A3. x = cotrans(x, -(90° - φ))  (to horizontal)
A4. Save horizon_lon; set x[0] = 0
A5. x = cotrans(x, -90°)  (to house circle)

For each cusp i (0-5):
  B0. x_cusp = [30°·i, 0, 1]
  B1. x_cusp = cotrans(x_cusp, 90°)
  B2. x_cusp[0] += horizon_lon
  B3. x_cusp = cotrans(x_cusp, 90° - φ)
  B4. x_cusp[0] += (ARMC - 90°)
  B5. λ = atan(tan(RA)/cos(ε)) with quadrant adjustment
```

---

## Natural Gradient (N)

**Historical Origin:** Modern simplified system.

**Conceptual Basis:** Equal houses starting from 0° Aries.

### Mathematical Formula

```
λᵢ = (i-1) × 30°
where i = 1, 2, ..., 12
```

House 1 = 0° Aries, House 2 = 0° Taurus, etc.

---

## APC (Y)

**Historical Origin:** Modern system (Ascendant-Parallel Circle).

**Conceptual Basis:** Based on the great circle parallel to the horizon passing through the Ascendant.

### Algorithm

Uses the declination of the Ascendant to define the house circle.

### Mathematical Formulas

**Ascensional difference of Ascendant:**
```
kv = arctan(tan(φ)·tan(ε)·cos(ARMC) / (1 + tan(φ)·tan(ε)·sin(ARMC)))
```

**Declination of Ascendant:**
```
δ_Asc = arctan(sin(kv) / tan(φ))
```

**House cusp calculation (for house n):**
```
For n < 8 (below horizon):
  k = n - 1
  a = kv + ARMC + π/2 + k·(π/2 - kv)/3

For n >= 8 (above horizon):
  k = n - 13
  a = kv + ARMC + π/2 + k·(π/2 + kv)/3

λ = atan2(tan(δ_Asc)·tan(φ)·sin(ARMC) + sin(a),
          cos(ε)·(tan(δ_Asc)·tan(φ)·cos(ARMC) + cos(a))
          + sin(ε)·tan(φ)·sin(ARMC - a))
```

---

## Summary Comparison

| System | Calculation Type | Latitude Independent | Works at Poles |
|--------|-----------------|---------------------|----------------|
| Placidus | Iterative | No | No |
| Koch | Iterative | No | No |
| Regiomontanus | Direct | No | Yes* |
| Campanus | Direct | No | Yes* |
| Porphyry | Direct | No | Yes |
| Equal | Simple | Yes | Yes |
| Whole Sign | Simple | Yes | Yes |
| Alcabitius | Direct | No | Yes* |
| Topocentric | Direct | No | Yes* |
| Morinus | Direct | Yes | Yes |
| Meridian | Direct | Yes | Yes |
| Vehlow | Simple | Yes | Yes |
| Carter | Direct | No | Yes |
| Gauquelin | Time-based | No | No |
| Horizontal | Complex | No | Yes* |
| Krusinski | Complex | No | Yes* |
| Natural | Simple | Yes | Yes |
| APC | Complex | No | Yes* |

*May require special handling near poles.

---

## References

1. Meeus, Jean. "Astronomical Algorithms" 2nd Edition (1998), Chapter 13
2. Swiss Ephemeris Documentation and Source Code
3. Hand, Robert. "Horoscope Symbols" (1981)
4. Holden, Ralph William. "The Elements of House Division" (1977)
5. IERS Conventions 2003 (IAU nutation models)
6. Polich, Wendel & Page, A.P. Nelson. "The Topocentric System of Houses" (1961)
7. Krusinski, Bogdan. Technical papers on Krusinski-Pisa system
