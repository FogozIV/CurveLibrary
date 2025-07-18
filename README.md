# CurveLibrary

CurveLibrary is a lightweight C++ library for generating and evaluating geometric curves, with a primary focus on clothoids (also known as Euler spirals). These curves are especially useful for path planning in robotics, offering continuous curvature transitions ideal for smooth trajectory generation.

## Purpose

This library is designed to:
- Provide efficient, portable clothoid curve generation
- Be compatible across platforms (Linux, Windows, microcontrollers like Teensy)
- Serve as a backend for trajectory previsualization and real-time embedded execution

## Origins and Acknowledgments

This project includes an adapted version of the original [Clothoids library](https://github.com/ebertolazzi/Clothoids) by:

- **Enrico Bertolazzi**  
- **Marco Frego**

Their excellent work on clothoid interpolation algorithms forms the mathematical and algorithmic foundation of this project. The core G² clothoid implementation is derived from their research publications:

- Bertolazzi, E., & Frego, M. (2013). *Fast and accurate G¹ fitting of clothoid curves*. arXiv:1305.6644.  
- Bertolazzi, E., & Frego, M. (2018). *On the G² Hermite Interpolation Problem with Clothoids*. Journal of Computational and Applied Mathematics, 341, 99–116. https://doi.org/10.1016/j.cam.2018.03.029

## Modifications

This version has been adapted with the following goals:
- Reduce dependencies and resource requirements for embedded systems
- Preserve numerical accuracy while improving portability
- Provide deterministic cross-platform behavior for path planning
- Allow integration with visualization tools (e.g., using `matplotlib`)

## License

This project uses the [BSD 2-Clause License](LICENSE), as required by the original Clothoids library. See the LICENSE file for details.

## Status

Currently supports G² clothoid interpolation with a simple API and builds on:
- Windows
- Linux
- Embedded platforms (tested on Teensy)

Further extensions (e.g., support for other curve primitives) are planned.
