**Plot boundary gradient vector**

The code here needs to be used in conjunction with the relevant methods in paper *Unified PID Control Analysis for Time-Delay Plants*. 

All of these coed is written in MATLAB. Among them, the code in program **Arg** is used to calculate the phase angle related to time delay, and continuity has already been taken into account. The Delta in program **plotGradientVector** and **plot_curves** can be directly obtained using the relevant methods in *Unified PID Control Analysis for Time-Delay Plants*. If necessary, it can also be obtained through program **computeDelta**, but it must be ensured that Delta has the same quantity as the gradient vector to guarantee numerical matching. Note that the function handles in **plot_curves** and **plotGradientVector** are the parametric equations of the relevant parametric curves or surfaces, while in **computeDelta**, *Fr* and *Fi* represent the real and imaginary parts of the characteristic equation when *s = iω*.
