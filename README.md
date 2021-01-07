## NVT SLLOD 

This code implements NVT-SLLOD for a system of soft-spheres in 2D. 
The temperature is controlled in the simulation using Nose-Hoover 
thermostat. 

To compile the code use 

`g++ main_nvt_sllod.cpp soft_sphere_nvt_sllod.cpp -o nvt.o`

Please use gcc version 7.3.0 or above

To run the code 

use 

`nvt.o [INPUT_FILE] [BOXh]`

![image](/PRESS_TEMP_NVT_SLLOD.png)
