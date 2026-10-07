# Sources and units

- Zhu, W., Allison, K. L., Dunham, E. M., and Yang, Y. (2020), Nature
  Communications 11, 4833, DOI 10.1038/s41467-020-18598-z. `paper.pdf`,
  `supplement.pdf`, and `published_fig*.png` are unmodified reference copies
  of CC BY 4.0 material. The project's generated figures are in `../figures/`.
- OSF DOI 10.17605/OSF.IO/9YGRP. `osf_manifest.json` records original names,
  sizes, SHA-256 values, and download URLs. `scripts/download_data.py` verifies
  all ten files. The nine `.mat` files are HDF5 and stay outside Git in `data/`.
  `data_information.pdf` is the authors' guide to figure/file correspondence.
- Scycle branch `Zhu_et_al_2020`, commit
  `d65e248900e193add0cc901ee7ee4d8ee7a2d50a`. `original_inputs/` contains the
  four unmodified case inputs. `SCYCLE_LICENSE` accompanies these MIT licensed
  files. The full source checkout is an ignored local reference, not part of
  the new solver. No original solver is silently invoked by the reproduction.
- `AGENTS.md` and the academic style guide were copied from the main branch of
  `bradlipovsky/vdv-damage`, observed commit
  `3d35f81420a744ca99360a6774d3bcbd4e123405`.
- Sources were accessed 2026-10-07. The OSF manifest is the authoritative record
  of archived file hashes; paper inputs and archived outputs are not assumed to
  be identical in undocumented details.

| Archived quantity | Stored units | Plot/solver conversion |
|---|---|---|
| depth | km | multiply by 1000 for C++ |
| time | s | divide by 31536000 for years |
| slip | m | unchanged |
| slip velocity | m/s | unchanged |
| effective normal stress | MPa | multiply by 1e6 for C++ |
| pore pressure | MPa | multiply by 1e6 for C++ |
| permeability | km² | multiply by 1e6 for m² |
| fluid flux | m/s | unchanged |

Important source distinctions:

- The archived effective stress includes the 1 MPa friction floor. The raw
  pressure-dependent permeability law does not use this floor.
- The supplied depth grid ends at 578.9526810302183 km, despite the declared
  `Lz = 500` in the input files. All archive files use this same grid.
- The T=1e7 s input uses flux 3.3764e-10 m/s, rather than the rounded 3.3e-10
  m/s in the paper caption.
- The T=1e10 s archive implies a deep reference permeability near 9.09e-16 m²;
  its input file lists 9.091e-16 m². Archive-assisted initial profiles preserve
  the recovered values, and evolve with the original input coefficients.
- The archive contains neither shear traction nor the aging-law state. Their
  reconstruction is approximate and recorded in `inputs/*_restart.json`.

The plotting scripts are new. Changed layouts, front-extraction conventions,
contour sampling, and time windows are documented in the report. Original
published images are never substituted for newly calculated or replotted panels.
