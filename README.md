# Educational_reverse-time-migration
Educational-level 2D acoustic reverse-time migration (RTM) code

My code works as a tool for seismologists to quickly understand how RTM imaging method works in synthetic dataset. This code provides a few options to test your imaging results:

(1) back-propagation solver (reconstruction vs. adjoint)

(Note that the classic RTM uses reconstruction to compute the receiver wavefield, which solves a boundary value problem. Nowadays, if we use the adjoint wave equation to get the receiver wavefield, I prefer to call it "adjoint migration" or "adjoint imaging", instead of RTM.)

(2) velocity model (true vs. smoothed)

(3) diagonal Hessian/source illumination (with vs. without)

Note that this code only gives you a flavor about RTM fundamental. Be careful if you use this code for your formal publications.

Contact chenglongduan.nju@gmail.com for any questions.
