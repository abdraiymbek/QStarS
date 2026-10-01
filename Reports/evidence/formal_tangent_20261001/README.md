# Report 19 transfer-tangent evidence

`nd20` and `nd39` contain the complete seven-family/six-direction audit,
including each input, native probe stdout and compressed numerical response.
`audit.json` records the exact source/library/compiler hashes and all step-size
results. `manifest.json` hashes every archived audit and numerical run file.

The selected gate sizes are explicit in the driver and identical on both grids.
Larger perturbations that cross face-interpolation limiter branches, and smaller
ones dominated by roundoff, remain in these records. A component pass does not
certify the native population Jacobian or a released cool solution.

`pytest.log` records the final regression suite. Initial diagnostic attempts,
including weak-source and oversized face-flux probes, remain in ignored
`results/formal_tangent`; their failure mechanisms and reproduced step-size
series are explained in Report 19.
