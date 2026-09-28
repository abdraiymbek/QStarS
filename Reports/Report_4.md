# Report 4 — Retraction of an audit interpretation

Date: 2026-09-28 UTC.

An initial audit note incorrectly asserted that Hassan et al.'s trigonometric
Eq 6 differed by a factor of two from its Appendix A cubic. Direct numerical
substitution disproves that assertion. At the critical mass, the cosine is
cos(pi/3)=1/2, so the displayed prefactor gives exactly the positive double
root. The small-mass limit also correctly tends to the BH-only Bondi radius.

The source registry and scientific documentation now state that the two
expressions agree. A regression compares the implemented cubic root against
the printed trigonometric expression as well as substituting into the cubic.
The implementation did not need a physical correction; the interpretation did.
The original erroneous note remains recoverable in Git commit f720887.

This correction does not affect the separately verified Ball/QUEST convention
difference or the Coughlin–Begelman entropy-minimum/printed-cutoff distinction.
