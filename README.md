# Educational_reverse-time-migration
Educational-level 2D acoustic reverse-time migration (RTM) code

My code works as a tool for seismologists to quickly understand how RTM imaging method works in synthetic dataset. This code provides a few options to test your imaging results:

(1) back-propagation solver (reconstruction vs. adjoint)

(Note that the classic RTM uses reconstruction to compute the receiver wavefield, which time-reversely 'recovers' the previous time step using the current time step wavefield (a boundary value problem). This is where the name RTM comes from.
However, if we use the adjoint wave equation to get the receiver wavefield, we compute the wavefield normally from 0 to T. The time-reverse concept in this case is indicated from the adjoint source data. Thus, I prefer to call it "adjoint migration" or "adjoint imaging" instead of RTM.)

(2) velocity model (true vs. smoothed)

(3) diagonal Hessian/source illumination (with vs. without)

The acoustic wave equation is a second-time-order PDE. The absorbing boundary uses exponential damping method. This code can give you a flavor about RTM fundamental. More detailed implementations (e.g., CPML boundary, first-order velocity-pressure scheme, HPC parallelism) are out-of-scope of this code.

Contact chenglongduan.nju@gmail.com for any questions.
