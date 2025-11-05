# 1 Overview

Multi-system and Multi-frequency Receiver Code Bias Analysis Software (MRCB) is an open-source GNSS data processing package specifically designed for extracting time-varying receiver code biases (RCBs). MRCB is developed in C/C++, which can be easily ported to different operating systems, such as Windows and Linux. It is a post-processing software, which can process multi-frequency data from GPS, Galileo, BDS, GLONASS, and QZSS. Using a least-squares filter (LSF) based on undifferenced and uncombined observations to estimate time-varying RCBs. Furthermore, it requires observation files and broadcast ephemeris, offering a simple operation process. Processing strategies for different data can be configured in the formatted configuration file. The main features of MRCB are as follows:

- Undifferenced and uncombined observations
- Support for both CDMA and FDMA
- Extract time-varying RCBs
- Multi-frequency (dual-frequency and above) data processing
- Least square filter
# 2 The project directory structure is as follows.
<img width="904" height="630" alt="filestructure" src="https://github.com/user-attachments/assets/698f766d-efd4-48da-a254-beb176a85334" />

# 3 Compile

## 3.1 Windows

Ensure that the complete compilation tool and CMAKE (version 3.10 or higher) is installed.

- Extract the files downloaded from GitHub
<img width="351" height="205" alt="image" src="https://github.com/user-attachments/assets/def70853-7f64-47cc-a916-5ff90b56b2a5" />


- Search "x64 Native Tools Command Prompt for VS 2022"
<img width="310" height="195" alt="image" src="https://github.com/user-attachments/assets/20068885-672c-4b0a-9a17-c2c9bfae0dc5" />

- Navigate to Project Directory and Compile. Execute these commands in sequence

① cd your_project_path ②mkdir build ③cmake -B build ④cmake --build build

## 3.2 Linux
Ensure that the complete compilation tool and CMAKE (version 3.10 or higher) is installed. Execute these commands in sequence: ① cd your_project_path ② mkdir build ③ cd build ④ cmake .. ⑤ make
<img width="799" height="566" alt="image" src="https://github.com/user-attachments/assets/6246aaa9-08d1-4c9c-9b16-1874f03f2c8a" />

