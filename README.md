# HoloRay

HoloRay is a C++ implementation of a **subpixel boundary refinement and area
estimation** module for closed contours on bitmap images.

The current version is implemented and tested with **Visual Studio 2015**
and uses a third-party KD-tree library for neighborhood queries.  
Porting to other compilers / environments should be possible by adjusting
project settings, but has not been tested.

---

## 1. Repository layout

HoloRay uses the following project structure:

    HoloRay.sln            # Visual Studio 2015 solution

    HoloRay/               # Main project
      data/                # Test data (images / boundaries / configs, etc.)
      KdTree/              # Third-party KD-tree implementation
      Parser/              # Simple parsers and small utilities
      Algo.h               # High-level algorithm entry points
      DataAnalysis.cpp
      DataAnalysis.h       # Analysis / statistics helpers
      FileProc.h           # File I/O helpers
      HoloRay.cpp          # Main driver / entry point
      ImgProc.cpp
      ImgProc.h            # Image / gradient processing helpers
      stdafx.cpp
      stdafx.h             # Precompiled header (VS2015)
      targetver.h          # Platform setup

The `data/` folder contains small test cases used for debugging and
experiments.

---

## 2. Requirements and build

### Environment

- IDE / compiler: **Microsoft Visual Studio 2015** (MSVC 14.x)  
- Language: C++ (uses basic C++11 features)  
- OS: Windows (tested on Windows 10)

### Build steps

1. Clone the repository:

       git clone https://github.com/<your-account>/<your-repo-name>.git

2. Open `HoloRay.sln` in **Visual Studio 2015**.

3. Select configuration:

   - `Release` for experiments / timing  
   - `Debug` for development and inspection  

4. Build the solution:

   - Menu: `Build → Build Solution`

If compilation fails, please check:

- Include path for `HoloRay/KdTree`  
- Runtime library settings (`/MT` vs `/MD`)  
- C++ language standard settings (C++11 or later)

---

## 3. Running HoloRay

1. Prepare your inputs in `HoloRay/data/`:

   - Test image(s)  
   - Closed pixel-level boundary  
   - Any additional files needed by your configuration  
   - (Optional) precomputed auxiliary data (e.g. charts / maps)

2. Adjust file paths and parameters in:

   - `HoloRay.cpp`  
   - `Algo.h`  
   - or related headers  

   Typical settings include input file names, output file names,
   and basic numerical parameters.

3. Build and run the `HoloRay` executable from Visual Studio
   (or from the build output directory).

4. Output (depending on configuration) may include:

   - Refined subpixel boundary samples  
   - Optional per-ray statistics / area estimates  

   Exact formats and file names can be checked in
   `HoloRay.cpp` and `DataAnalysis.cpp`.

---

## 4. Third-party code

The `HoloRay/KdTree` directory contains a KD-tree implementation
written by an external author. Please refer to the license / README
inside `KdTree/` for details and licensing terms. It is used only as
a search acceleration structure.

All other source files in this repository are covered by the license
specified in the top-level `LICENSE` file (to be chosen by the
repository owner, e.g. MIT / BSD / Apache-2.0).

---

## 5. Citation and usage

This code is a **research prototype**:

- Only Visual Studio 2015 on Windows has been tested end-to-end.  
- Interfaces and file formats may change as the project evolves.  
- The focus is on clarity of the subpixel refinement pipeline,
  not on production-grade optimization.

Bug reports and pull requests are welcome.
