# FIB-SEM Toolkit

A Fiji/ImageJ plugin suite for quantitative analysis of Focused Ion Beam - Scanning Electron Microscopy (FIB-SEM) images.

## Features

- **Read .dat Files** — Parse custom binary FIB-SEM data files with full
  header metadata, spatial/temporal calibration, and multi-channel support.
- **Noise Statistics Analysis** — Characterize image noise via variance–mean
  regression, estimating system gain and signal-to-noise ratios.
- **Contrast Analysis** — Quantify contrast between material phases using
  CDF-based peak identification on gradient-filtered intensity distributions.
- **Edge Transition Analysis** — Measure edge sharpness by detecting gradient
  peaks, extracting intensity profiles along gradient directions, and fitting
  a directional ellipse to characterize resolution anisotropy.
- **Gradient Mapping** — Compute gradient magnitude maps with optional
  smoothing and intensity normalization.
- **Threshold Calculation** — Estimate intensity thresholds from histogram
  quantiles with PDF/CDF visualization.

## Requirements

- [Fiji](https://fiji.sc/) (ImageJ2 distribution)
- Java 8+
- Maven (for building from source)

## Installation

### From Source

1. Clone the repository
2. Build with Maven:
mvn clean install


3. Copy the resulting JAR from `target/` into your Fiji `plugins/` directory.
4. Restart Fiji.

### Manual

Copy a pre-built `fib_sem-<version>.jar` into your Fiji `plugins/` directory
and restart Fiji.

## Usage

All plugins are available under **Plugins > FIB-SEM** in the Fiji menu bar.

### Recommended Workflow

1. **Read .dat** — Open a FIB-SEM `.dat` file. Metadata is attached via
`Image > Show Info`.
2. **Noise Statistics** — Run on the opened image to estimate system gain,
dark count (I0), and SNR. Results are persisted on the image for
downstream use.
3. **Contrast Analysis** — Computes contrast using the dark count from step 2
(or a manually entered value).
4. **Edge Transition Analysis** — Measures edge sharpness and directional
anisotropy.
5. **Gradient Map / Thresholds** — Auxiliary tools for visualization and
threshold estimation.

### ROI Support

All analysis plugins respect the active ROI. Draw a rectangular ROI before
running a plugin to restrict analysis to a specific region.

### Parameter Persistence

Plugin parameters are stored on each image and persist across invocations
within the same Fiji session. This allows iterative refinement without
re-entering values.

### Figure Export

Most plugins include a "Save as figure" option that writes titled PNG files
to the source image's directory.

## Project Structure
```bash
src/main/java/com/shtengel/fib_sem/
├── core/           # Analysis algorithms (stateless, static methods)
│   ├── ContrastAnalyzer.java
│   ├── DatFileProcessor.java
│   ├── EdgeTransitionAnalyzer.java
│   ├── GradientMapAnalyzer.java
│   ├── NoiseStatisticsAnalyzer.java
│   └── ThresholdAnalyzer.java
├── data/           # Immutable result containers
│   ├── ContrastData.java
│   ├── EdgeTransitionData.java
│   ├── FileData.java
│   ├── GradientMapData.java
│   ├── HeaderData.java
│   ├── NoiseStatisticsData.java
│   └── ThresholdData.java
├── plugins/        # Fiji Command plugins (UI layer)
│   ├── Contrast.java
│   ├── EdgeTransitions.java
│   ├── GradientMap.java
│   ├── NoiseStatistics.java
│   ├── ReadDatFile.java
│   └── Thresholds.java
└── util/           # Shared utilities
├── EllipseFitter.java
├── FigBuilder.java
├── ImageResolver.java
└── ParamPersister.java
```


## Architecture

The codebase follows a three-layer pattern:

- **Plugins** handle UI dialogs, parameter persistence, visualization, and
  figure export. They delegate all computation to core analyzers.
- **Core analyzers** are stateless classes with static methods that accept
  ImageJ `ImageProcessor` inputs and return typed data containers.
- **Data containers** are immutable records holding analysis results.

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.