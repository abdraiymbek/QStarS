# MESA and production fidelity

MESA was absent from the initial environment. Its documented installation
requires a version-matched SDK and about 20 GB of disk space. This Mac has
enough storage and memory, but no MESA build or stellar model was executed.
These are outstanding tasks, not evidence that MESA cannot run here.

MESA-QUEST was inspected at commit
`671eb9b898b28c0d02a6a2455dacc0f5b269e280`. Its README associates the early
Bondi model with MESA r24.03.1, the saturated/wind model with r24.08.1 or
r25.12.1, and its colors model with r25.12.1. It is cached, not vendored. No
upstream executable was run.

Full scientific evolution should use established MESA microphysics with
independently audited QStarS boundary accounting. Remaining work includes:

1. Pin and build MESA/SDK/QUEST and reproduce an unmodified stellar benchmark.
2. Audit startup, retries, restart storage, and core relaxation for BH/cavity
   variables. Upstream module-level BH state motivates this audit; it is not
   a demonstrated upstream retry bug.
3. Add independent boundary enthalpy/entropy streams and matching potential work.
4. Establish conservative restart/remesh/retry tests for that actual backend.
5. Reproduce multiple full evolution and advective-wind benchmarks.

No placeholder MESA adapter is advertised as implemented.
[Installation](https://docs.mesastar.org/en/latest/installation.html),
[energy equations](https://docs.mesastar.org/en/latest/reference/controls.html#energy-eqn-option),
[MESA-QUEST](https://github.com/andysantarelli/MESA-QUEST).
