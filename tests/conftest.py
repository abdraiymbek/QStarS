import os
import sys

# The MESA SDK gfortran driver passes the host macOS version to the system linker/assembler;
# on hosts newer than the SDK toolchain this is rejected ("invalid version number in
# -mmacosx-version-min"). backends/cmfgen/build.sh pins the same target for native builds.
if sys.platform == 'darwin':
    os.environ.setdefault('MACOSX_DEPLOYMENT_TARGET', '14.0')
