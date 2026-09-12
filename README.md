# REACH South Asia - Interactive Web Model

🌍 **[Live Application](https://mrafiuddin50.github.io/REACH-wasm/)**

REACH (Reduced-form Evaluation of Air quality Control policies for Health) is a state-of-the-art atmospheric simulation model designed to track the dispersion and chemistry of air pollution (PM2.5) across South Asia. 

This repository contains the completely serverless, client-side web application version of the REACH model. It allows users to run complex atmospheric dispersion and chemical thermodynamic equilibrium scenarios instantaneously within their web browsers.

---

## 🔬 How REACH Works

The core REACH model simulates how primary emissions (like SO2, NOX, and NH3) from distinct sectors mix into the atmosphere and react over long distances to form deadly PM2.5. The model is composed of two primary engines:

### 1. Atmospheric Physics (Gaussian Plume Dispersion)
In the traditional Python/R model, REACH evaluates hundreds of millions of spatial relationships using a Gaussian Plume dispersion model driven by highly granular 2022 meteorological data (wind speed, wind direction, temperature, and planetary boundary layer heights). 

### 2. Secondary Aerosol Chemistry
Once primary pollutants disperse, REACH applies a rigorous thermodynamic equilibrium chemistry model (similar to ISORROPIA). It evaluates the complex interactions between Sulfates (SO4), Nitrates (NO3), and Ammonium (NH4) based on the availability of Ammonia (NH3) gas, while also tracking Secondary Organic Aerosols (SOA) and primary Black Carbon (BC) and Organic Carbon (OC).

---

## ⚡ The Web Architecture

To convert a massive, compute-heavy atmospheric model into an interactive application that runs at 60fps on a mobile phone without a server, we engineered a completely novel architecture.

### Pre-Compiled Source-Receptor Matrices
Instead of solving the complex Gaussian Plume equations in real-time, the meteorological physics are pre-compiled in Python. We condensed the dispersion data for all 1,136 South Asian districts into lightweight **Source-Receptor (SR) matrices**. These matrices act as a mathematical "lookup table" that encode precisely how one ton of emissions in District A will disperse into District B over the course of a year.

### The Raw Binary Engine (Bypassing WASM)
Initially, this application utilized WebAssembly (WASM) and Pyodide to run the original Python matrices in the browser. However, initializing the heavy Python WASM engine and decompressing the Numpy files caused significant loading delays.

To achieve truly instantaneous loading and execution, the Python `numpy` DataFrames were converted into flat, C-contiguous **Float32 Binary Arrays**. 
1. **Zero Initialization Overhead:** The browser fetches these raw `.bin` files and maps them directly into native memory using JavaScript `TypedArrays`, bypassing WASM/Pyodide entirely.
2. **Native Matrix Multiplication:** When a user adjusts an emission slider, the browser's native JavaScript V8 Engine executes the massive matrix multiplications in a highly optimized `O(N^2)` loop at native C++ speeds.
3. **Client-Side Chemistry:** The thermodynamic equilibrium chemistry simulation is calculated natively in JavaScript for every grid cell, outputting the final PM2.5 concentrations in milliseconds.

### Serverless & Private
Because the math is handled via binary arrays and native browser JIT compilation, **there is no backend server**. Everything happens securely, privately, and instantaneously on your device.
