# terSIDH: Constant terSIDH (p128 by default)

## Build

To build with default parameters, simply run:

    make

Configurable build arguments:

- BITS Selects a supported parameter set directory (default: 128 → `p128/`).
- UINT_IMPL Unsigned‑integer backend implementation.
- FP_IMPL Finite‑field backend implementation.

Examples:

    make BITS=1024

builds using the `p1024/` parameter set (if present), and:

    make UINT_IMPL=uint.c FP_IMPL=fp.c

uses the generic C arithmetic instead of assembly.

## Adding parameter sets

Create a directory `p${BITS}/` modelled after `p128/`, optionally with
specialized `uint*.s` and `fp*.s` for assembly arithmetic. Otherwise the
generic C implementations will be used.

## Targets

- `make` Builds the shared library (`libtersidh.so`) and the demo `main`.
- `make test` Builds a simple test program.
- `make clean` Removes build artifacts.

prints aggregate timing and operation counts for 1000 iterations.

## About this implementation

This codebase implements a constant‑time variant of the Montgomery isogeny
toolchain for terSIDH‑style flows. Key elements:

- Per‑path cofactor selection for kernel reduction (real/dummy suffixes) to ensure
  the kernel point has exact order `primes[i]` at every step (including degree‑4).
- Fixed‑iteration Montgomery ladder is used in the key‑independent “prep” phase to
  mitigate timing variance from differing scalar bitlengths.
- Constant‑time conditional moves (cmov) are used for path selection and state updates.

The provided `main` demonstrates key generation and prints a compact table across
several private‑key distributions. The demo can be adapted for integration or further
measurement as needed.

## References

This implementation makes use of ideas from these references:

- W. Castryck, T. Lange, C. Martindale, L. Panny, J. Renes:
  CSIDH: An efficient post‑quantum commutative group action.
  Asiacrypt 2018, https://ia.cr/2018/383

- M. Meyer, S. Reith:
  A faster way to the CSIDH.
  Indocrypt 2018, https://ia.cr/2018/782

- D. J. Bernstein, T. Lange, C. Martindale, L. Panny:
  Quantum circuits for the CSIDH: Optimizing quantum evaluation of isogenies.
  Eurocrypt 2019, https://ia.cr/2018/1059

- M. Meyer, F. Campos, S. Reith:
  On Lions and Elligators: An efficient constant‑time implementation of CSIDH.
  PQCrypto 2019, https://ia.cr/2018/1198

- D. J. Bernstein, B.-Y. Yang:
  Fast constant‑time gcd computation and modular inversion.
  CHES 2019, https://ia.cr/2019/266

- S. Kim:
  Performance and Efficiency Evaluation of M‑SIDH.
  IEEE Access, 2024. https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10786215
