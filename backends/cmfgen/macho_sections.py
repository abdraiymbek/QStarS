#!/usr/bin/env python3
"""Path-independent fingerprint of a native CMFGEN build (build 33 onwards).

The sha256 of `cmfgen_dev.exe` depends on where and when it was linked: ld writes a debug map of N_OSO stabs that
name each archive member by its absolute path (.../build_NN/cur_cmf/lib/libsubs.a(member.o)) with the member's
timestamp, plus LC_UUID and the ad-hoc code signature. None of these is inside a section. The sha256 of the file
content of every section (zero-fill sections: their size) and of the changed archive members is therefore the same
for the same sources, compiler and flags in any directory.

    macho_sections.py EXE [LIBSUBS MEMBER.o ...]   ->  JSON {sections: {...}, objects: {...}}
"""
import hashlib
import json
import struct
import subprocess
import sys
from pathlib import Path

LC_SEGMENT_64 = 0x19
ZEROFILL = {0x1, 0xc, 0x12}          # S_ZEROFILL, S_GB_ZEROFILL, S_THREAD_LOCAL_ZEROFILL


def section_sha256(path):
    """{'SEG,sect': sha256 of the section's file bytes, or 'zerofill:SIZE'} of a thin 64-bit Mach-O file."""
    b = Path(path).read_bytes()
    magic, _, _, _, ncmds = struct.unpack_from('<IiiII', b, 0)
    if magic != 0xfeedfacf:
        raise ValueError(f'{path}: not a thin 64-bit little-endian Mach-O file')
    out, off = {}, 32
    for _ in range(ncmds):
        cmd, size = struct.unpack_from('<II', b, off)
        if cmd == LC_SEGMENT_64:
            nsects = struct.unpack_from('<I', b, off + 64)[0]
            s = off + 72
            for _ in range(nsects):
                name = b[s:s + 16].rstrip(b'\0').decode()
                seg = b[s + 16:s + 32].rstrip(b'\0').decode()
                sz, foff = struct.unpack_from('<QI', b, s + 40)
                flags = struct.unpack_from('<I', b, s + 64)[0]
                key = f'{seg},{name}'
                if key in out:
                    raise ValueError(f'{path}: duplicate section {key}')
                out[key] = f'zerofill:{sz}' if (flags & 0xff) in ZEROFILL or foff == 0 else \
                    hashlib.sha256(b[foff:foff + sz]).hexdigest()
                s += 80
        off += size
    return out


def member_sha256(archive, member):
    """sha256 of one member of a static library, as stored (ar -p)."""
    r = subprocess.run(['ar', '-p', str(archive), member], capture_output=True, check=True)
    if not r.stdout:
        raise ValueError(f'{archive}: no member {member}')
    return hashlib.sha256(r.stdout).hexdigest()


if __name__ == '__main__':
    exe = sys.argv[1]
    rec = dict(sections=section_sha256(exe))
    if len(sys.argv) > 3:
        rec['objects'] = {m: member_sha256(sys.argv[2], m) for m in sys.argv[3:]}
    print(json.dumps(rec, indent=1))
