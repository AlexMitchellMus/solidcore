# solidcore

A B-rep solid modeling kernel in safe, deterministic, wasm-compatible Rust: a Rust port
of NVIDIA's SMLib (from [USD BRep](https://github.com/NVIDIA-Omniverse/usd-brep),
Apache-2.0). Not affiliated with or endorsed by NVIDIA.

**Status: 1.0.1.** The version tracks the upstream release it is ported from and
tested against, usd-brep v1.0.1. All of SMLib is ported:
- primitives, booleans, fillets and chamfers, offsets and shells, sweeps, healing, and
  mass properties;
- a unit is green only after differential testing against the C++ original confirms
  identical topology and values within 1e-9 relative, with any per-key absolute floors
  documented in its report.

The public API mirrors SMLib's and may still change between minor versions while
integrations settle.

Current unit results, blockers and primitive checks:
<https://tangency-labs.github.io/solidcore/>.

Progress, design and how to contribute: <https://github.com/tangency-labs/solidcore>
(see `plan/units/`, `plan/PROGRESS.md` and `AGENTS.md`).

Licensed under Apache-2.0; see `LICENSE`, `NOTICE` and `ATTRIBUTION.md`.
