---
name: Bug report
about: Report something that does not work as documented
labels: bug
---

**What happened**
A short description of the problem and the error message, if any.

**Minimal reproducible example**
Please include code that runs in a fresh R session, ideally with
`reprex::reprex()`. Use small dimensions, `dt_emulation`, `dt_calibration`, or
simulated data.

```r
library(spDBL)

```

**What you expected**

**Session information**
Paste the output of `sessionInfo()`. For calibration problems, also give the
rstan and StanHeaders versions and your C++ toolchain.
