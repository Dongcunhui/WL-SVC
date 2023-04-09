# Wavelet-Based Learned Scalable Video Coding (WL-SVC)
## Overall
This repo provides the official implementation of "Wavelet-Based Learned Scalable Video Coding".

The overall framework of WL-SVC is shown in the figure below. It consists of two parts: traditional modules (blue parts) and neural network modules (orange parts). Traditional modules come from [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) which is implemented in C, and neural network modules are implemented in python by us.
<img src="https://user-images.githubusercontent.com/48936648/150902487-f2288ab0-0a8d-4cb9-90b9-8b918dd59854.png" width="500px">

## How To Train
The code is implemented by pytorch and requires torch version >= 1.6.  

The residual coding module and inverse wavelet transform module are trainable, and the code is in train folder.

The input are original frames, motion vector and mask. You can obtain the motion vector and mask by modifying [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) yourself or refer to my modifications in the Interframe EZBC folder.

## How To Test
The test consists of two steps,
1. Use the modified Interframe EZBC-JP2K code to obtain the motion vector and the index of the reference frame, and use Interframe EZBC-JP2K to encode them into a code stream. (You can modify [Interframe EZBC-JP2K](https://ecse.rpi.edu/interframevideocoding/) yourself or refer to my modifications in the Interframe EZBC folder.)

2. Use the test.py code in the test folder to code a video. The input is the YUV component of the original video, the motion vector and the mask of the reference frame. The input data of an example is stored in the [network disk](https://drive.google.com/drive/folders/1cGloGAZZtUtqbWPC5-SZBsm8tD9RUm1U?usp=sharing), and after downloading, place it in the test folder. The trained model is stored in the [network disk](https://drive.google.com/drive/folders/1wVlfJ1tH1UdyttPOwYA2lURqYHJm5hQK?usp=share_link), including the entropy coding model and the post-processing model. After downloading, place them in the test folder. (The code is implemented by pytorch and requires torch version >= 1.6.)
