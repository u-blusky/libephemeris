# TODO - Improvements and Technical Debt

All items completed as of March 2026.

---

## Completed

- ~~#1~~ Thread safety — `_STATE_LOCK` RLock on critical setters/getters in state.py
- ~~#2~~ Bare `except Exception` — 101 -> 0 across entire codebase (all narrowed to specific types)
- ~~#3~~ MeeusRangeError — marked deprecated
- ~~#4~~ houses.py dedup — `_init_cardinal_cusps` (9x) + `_set_opposite_cusps` (6x) extracted
- ~~#5~~ constants.py `__all__` — 785 public names, `annotations` excluded
- ~~#6~~ PyPI Packaging — `base_core.leb` bundled in wheel, auto-discovered
- ~~#7~~ Download command — `download_leb2_for_tier()` in download.py, 12 files in DATA_FILES
- ~~#8~~ LEB1 regenerated — Pluto 64d/deg11, Uranians 256d/deg7
- ~~#9~~ GitHub Release data-v2 — 12 LEB2 files published
- ~~#10~~ Validation suite — 8/8 PASS, 61,547 checks, ALL PASS
- ~~#11a~~ Uranian geocentric — added geocentric path in `planets.py`
- ~~#11b~~ Sun heliocentric — returns (0,0,0) correctly now
- ~~#11c~~ True Node distance — full osculating orbit calc (2000x improvement)
- ~~#12~~ Skyfield TypeError — retry + radec removal, cache clear
- ~~#13~~ Developer CLI split — `leph` now has its own `status` command and a development bootstrap `download all` that fetches prerequisites and source inputs, while `libephemeris` stays focused on end-user runtime setup
