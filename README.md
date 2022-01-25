# Wavelet-Based-Learned-Scalable-Video-Coding (WL-SVC)
## Overall
This repo provides the official implementation of "Wavelet-Based Learned Scalable Video Coding". It is just a preliminary version, and details will be added later.  
The overall framework of WL-SVC is shown in the figure below. It consists of two parts: traditional modules (blue parts) and neural network modules (orange parts). Traditional modules come from [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) which is implemented in C, and neural network modules are implemented in python  by us.
<img src="https://user-images.githubusercontent.com/48936648/150902487-f2288ab0-0a8d-4cb9-90b9-8b918dd59854.png" width="500px">

## How To Train
1. The residual coding module and inverse wavelet transform module are trainable, and the code in . The input 
1. Use [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) to get metion vector and mask.
2. input original frame, MV and mask into python.
3. the MV bits is grenerate by [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) and the residual bits is generate by python.
4. Have fun!
## How To Test
