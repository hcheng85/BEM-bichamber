#!/usr/bin/env python3
"""Minimal read-only memory mapping for uncompressed MATLAB v5 numeric arrays.

The circular trajectory files store very large ``miMATRIX`` elements without
``miCOMPRESSED`` wrappers.  MATLAB uses Fortran order, so ``(component, fish,
time)`` frames are contiguous and can be processed sequentially without ever
materializing a complete variable.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np


MI_INT8 = 1
MI_UINT8 = 2
MI_INT16 = 3
MI_UINT16 = 4
MI_INT32 = 5
MI_UINT32 = 6
MI_SINGLE = 7
MI_DOUBLE = 9
MI_INT64 = 12
MI_UINT64 = 13
MI_MATRIX = 14
MI_COMPRESSED = 15

DTYPES = {
    MI_INT8: "i1", MI_UINT8: "u1", MI_INT16: "i2", MI_UINT16: "u2",
    MI_INT32: "i4", MI_UINT32: "u4", MI_SINGLE: "f4", MI_DOUBLE: "f8",
    MI_INT64: "i8", MI_UINT64: "u8",
}


def align8(n: int) -> int:
    return (n + 7) & ~7


@dataclass(frozen=True)
class NumericArrayInfo:
    path: Path
    name: str
    shape: tuple[int, ...]
    dtype: np.dtype
    data_offset: int
    data_nbytes: int

    def memmap(self) -> np.memmap:
        return np.memmap(
            self.path, dtype=self.dtype, mode="r", offset=self.data_offset,
            shape=self.shape, order="F",
        )


def _tag(handle, endian: str) -> tuple[int, int, int, int]:
    """Return ``(type, nbytes, data_offset, next_offset)`` for either tag form."""
    pos = handle.tell()
    raw = handle.read(8)
    if len(raw) != 8:
        raise EOFError
    first, second = struct.unpack(endian + "II", raw)
    if first >> 16:
        dtype = first & 0xFFFF
        nbytes = first >> 16
        if nbytes > 4:
            raise ValueError(f"Invalid small-data element of {nbytes} bytes")
        return dtype, nbytes, pos + 4, pos + 8
    return first, second, pos + 8, pos + 8 + align8(second)


def inspect_numeric_arrays(path: str | Path) -> dict[str, NumericArrayInfo]:
    path = Path(path)
    result: dict[str, NumericArrayInfo] = {}
    with path.open("rb") as handle:
        header = handle.read(128)
        if not header.startswith(b"MATLAB 5.0 MAT-file"):
            raise ValueError(f"Not a MATLAB v5 file: {path}")
        indicator = header[126:128]
        endian = "<" if indicator == b"IM" else ">" if indicator == b"MI" else None
        if endian is None:
            raise ValueError(f"Unknown endian indicator {indicator!r}")
        endian_prefix = "<" if endian == "<" else ">"

        file_size = path.stat().st_size
        while handle.tell() + 8 <= file_size:
            outer_pos = handle.tell()
            try:
                outer_type, outer_nbytes, outer_data, next_outer = _tag(handle, endian)
            except EOFError:
                break
            if outer_type == MI_COMPRESSED:
                raise ValueError(
                    f"{path} contains miCOMPRESSED data; it cannot be memory mapped"
                )
            if outer_type != MI_MATRIX:
                handle.seek(next_outer)
                continue

            handle.seek(outer_data)
            # Array flags.
            typ, nbytes, off, next_element = _tag(handle, endian)
            handle.seek(next_element)
            # Dimensions.
            typ, nbytes, off, next_element = _tag(handle, endian)
            if typ not in (MI_INT32, MI_UINT32):
                raise ValueError("Unexpected MATLAB dimensions element")
            handle.seek(off)
            dims = struct.unpack(endian + "I" * (nbytes // 4), handle.read(nbytes))
            handle.seek(next_element)
            # Name (often stored with MATLAB's compact small-data tag).
            typ, nbytes, off, next_element = _tag(handle, endian)
            handle.seek(off)
            name = handle.read(nbytes).decode("latin1")
            handle.seek(next_element)
            # Numeric payload. Structs/cells are skipped.
            typ, nbytes, off, next_element = _tag(handle, endian)
            expected = int(np.prod(dims, dtype=np.int64))
            if typ in DTYPES and nbytes == expected * np.dtype(DTYPES[typ]).itemsize:
                dtype = np.dtype(endian_prefix + DTYPES[typ])
                result[name] = NumericArrayInfo(path, name, tuple(dims), dtype, off, nbytes)
            handle.seek(next_outer)
    return result


def open_state_rdot(path: str | Path) -> tuple[np.memmap, np.memmap, dict[str, NumericArrayInfo]]:
    info = inspect_numeric_arrays(path)
    missing = {"state", "rdot"} - set(info)
    if missing:
        raise KeyError(f"Missing numeric arrays {sorted(missing)} in {path}")
    return info["state"].memmap(), info["rdot"].memmap(), info
