# SFND_Radar_Target_Generation_and_Detection

This project is part of the sensor fusion nanodegree. Here a RADAR transmission is simulated for a given range and target velocity of the vehicle being tracked.
Subsequently the reflection of that transmission is processed and the distance and velocity of the target object is calculated and plotted.

# 2D CFAR

CFAR stands for constant false alarm rate. The purpose of CFAR is to reliably identify peaks in a surface while maintaining a predictable rate of false alarms.
It works by setting the threshold for peak detection to be an offset from the average values in the vicinity of the cell being evaluated for whether it is a peak or not.
The value of this offset threshold is determined experimentally for the purpose at hand. The number of cells to evaluate in the vicinity of the target cell is also an empirical decision.
The way the averaging is done is also an implementation decision. In the current project a mean of all training cells around the target cell, separated a by a set of guard cells is taken.
The size of the training grid is 11x11. A 5x5 grid within this 11x11 grid is a set of guard cells. This gives a total of 96 training cells for each target cell.
The offset threshold value used is 6dB. This implies that the target cell needs to be 4 times as powerful as any surrounding cells to be picked as a peak.
In the experiments with the above values, reliable velocity and range measurements were seen.

# Selection of Training and Guard cells and offset threshold

For the current problem, various training and guard cell sizes were tried, at lower values of training cells, we see that the peak of the distribution can get suppressed as the gradient at the
peak is lower than the gradient leading up to the peak, implying that there is a region of relatively flat response near the peak, which will suppress the peak if the training cells are too close
to the target cell. Additionally, having a low threshold tends to pick up more false positives. Having cycled through various values, a training grid (including the target cell and guard cells) of 11x11
and a guard grid (including target cell) of 5x5 appears to be ideal along with a offset threshold of 6dB.

# Steps to suppress non-thresholded cells

Due to the nature of the grid, the cells at the periphery of the grid cannot be evaluated in the same way as all the other cells.
There are a couple of options available:

1. Completely ignore those cells and just assume that there is no peak there. This will cause issues with detection of range and velocities near the extremes
2. Use different grid sizes closer to edges to get some results for those regions, but assign them a confidence score to reflect the different parameters used.

In this project option 1 was used.
