# Automotive Cruise Control Simulator (PID)

An advanced discrete-time Cruise Control simulator designed for a 1500 kg vehicle. This project bridges theoretical automatic control systems with real-world embedded software architecture.

## 🛠️ Tech Stack & Architecture
* **Core Simulation (Embedded Layer):** Pure C implementation of a discrete-time PID controller.
* **Data Analysis Pipeline:** Python (Pandas & Matplotlib) to parse telemetria logs and plot performance metrics.

## 🧠 Key Control Engineering Features
* **Anti-Windup (Back-Calculation):** Prevents integral saturation due to physical actuator limits (Engine max force limited to 5000 N).
* **Slew-Rate Limiter:** Smooths the reference input signal (ramp profile) to prevent mechanical shock and simulate realistic acceleration.
* **Performance Evaluation:** Automated calculation of **ISE** (Integral Squared Error), **IAE** (Integral Absolute Error), Rise Time, and Max Overshoot.

## 🚀 How to Run the Project
1. Compile and run the C simulation to generate the dataset:
   ```bash
   gcc cruise_control.c -o cruise_control.exe
   .\cruise_control.exe
