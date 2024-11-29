"""


---

"""

from math import pi

import numpy as np
import math
from matplotlib import pyplot as plt
from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

""" **USER INPUTS**

These parameters are typically on the user interface of the scanner computer console 
"""
xyanlge=0
xy_ang = round(xyanlge * pi / 180, 3)
nsa = 1  # Number of averages
n_slices = 1  # Number of slices
Nx = 256
Ny = 256
fovx = 220.e-3  # mm
fovy = 220.e-3
slice_thickness = 2.e-3  # s
slice_gap = 0.e-3  # s
rf_flip = 90  # degrees
rf_offset = 0
apodization=1
print('User inputs setup')

"""
Set the hardware limits and initialize sequence object
"""

system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', 
              rf_ringdown_time=10e-6, grad_raster_time=10.e-6,
              rf_dead_time=100e-6)
seq = Sequence(system)

""" TIME CONSTANTS """

TE = 7.1e-3  # s
TR = 4.5  # s
TI= 7.e-3
readout_time = 2.5e-3
pre_time = 8e-4  # s

""" RF """

flip90 = round(rf_flip * pi / 180, 3) #round to three decimal digits
flip180 = 180 * pi / 180
rf90, gz90, _ = make_sinc_pulse(flip_angle=flip90, system=system, duration=1.99e-3, ##Calculate amplitude from flip=360 gamma \int_{t}^{t+duration}B1(t) dt
                                slice_thickness=slice_thickness, apodization=apodization, 
                                return_gz=True,time_bw_product=2, use = 'excitation')
rf180, gz180, _ = make_sinc_pulse(flip_angle=flip180, system=system, return_gz=True,
                                  duration=1.89e-3, 
                                  slice_thickness=slice_thickness, 
                                  apodization=apodization, 
                                time_bw_product=2, use = 'refocusing')
rf_inv,gz_inv,_ = make_sinc_pulse(flip_angle=flip180, system=system, return_gz=True,
                                  duration=2.51e-3, 
                                  slice_thickness=slice_thickness, 
                                  apodization=apodization, 
                                time_bw_product=2, use = 'inversion')

"""
Readout gradients and related events
"""
delta_kx = 1 / fovx
delta_ky = 1 / fovy
k_width_x = Nx * delta_kx
k_width_y = Ny * delta_ky
gx = make_trapezoid(channel='x', system=system, flat_area=k_width_x, 
                    flat_time=readout_time)
adc = make_adc(num_samples=Nx, duration=gx.flat_time, delay=gx.rise_time)
print(gx)

""" REPHASE AND REPHASE """

phase_areas = (np.arange(Ny) - (Ny / 2)) * delta_ky
gz_reph = make_trapezoid(channel='z', system=system, area=-gz90.area / 2,
                         duration=readout_time / 2)
gx_pre = make_trapezoid(channel='x', system=system, flat_area=k_width_x / 2, 
                        flat_time=readout_time / 2)
gy_pre = make_trapezoid(channel='y', system=system, area=phase_areas[-1], 
                        duration=readout_time / 2)



""" 
Echo time (TE) and repetition time (TR). Here, TE is broken down into `delay1` and `delay2`.
"""
delay_TI = TI - calc_duration(rf_inv) / 2 - calc_duration(rf90) / 2

delay_TE1=TE/2 - calc_duration(rf90) / 2 -calc_duration(gx_pre) - calc_duration(rf180) / 2

delay_TE2=TE/2 - calc_duration(rf180) / 2 - calc_duration(gx)/2 
 
delay_TR = TR -TI - TE - calc_duration(gx)/2


delay_TR = make_delay(delay_TR)
delay_TE1 = make_delay(delay_TE1)
delay_TE2 = make_delay(delay_TE2)
delay_TI = make_delay(delay_TI)
print(f'delay_TI: {delay_TI}')
print(f'delay_TE: {delay_TE1}')
print(f'delay_TE: {delay_TE2}')
print(f'delay_TR: {delay_TR}')

"""
Construct sequence for one phase encode and multiple slices
"""

# Prepare RF offsets. This is required for multi-slice acquisition
delta_z = n_slices * slice_gap
z = np.linspace((-delta_z / 2), (delta_z / 2), n_slices) + rf_offset
#print(phase_areas)
for k in range(nsa):  # Averages
  for j in range(n_slices):  # Slices
    # Apply RF offsets
    freq_offset = gz90.amplitude * z[j]
    rf90.freq_offset = freq_offset

    freq_offset = gz180.amplitude * z[j]
    rf180.freq_offset = freq_offset

    for i in range(Ny):  # Phase encodes
        seq.add_block(rf_inv,gz_inv)
        seq.add_block(delay_TI)
        seq.add_block(rf90, gz90)
        gy_pre = make_trapezoid(channel='y', system=system, 
                                area=phase_areas[i], duration=readout_time / 2)
        seq.add_block(gx_pre, gy_pre, gz_reph)
        seq.add_block(delay_TE1)
        seq.add_block(rf180, gz180)
        seq.add_block(delay_TE2)
        seq.add_block(gx, adc)

        seq.add_block(delay_TR)

""" PLOTTING TIMING DIAGRAM """

ok, error_report = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print('Timing check passed successfully')
else:
    print('Timing check failed. Error listing follows:')
    [print(e) for e in error_report]
seq.plot(time_range=(TR-0.05, TR+0.2))

""" GENERATING `.SEQ` FILE 
Uncomment the code in the cell below to generate a `.seq` file and download locally.
"""


#seq.plot()



# Trajectory calculation and plotting
#ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()
#time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
#plt.figure(figsize=[12, 7])
#plt.plot(time_axis, ktraj.T)  # Plot the entire k-space trajectory
#plt.plot(t_adc, ktraj_adc[0], '.')  # Plot sampling points on the kx-axis
#plt.figure(figsize=[12, 7])
#plt.plot(ktraj[0], ktraj[1], 'b')  # 2D plot
#plt.axis('equal')  # Enforce aspect ratio for the correct trajectory display
#plt.plot(ktraj_adc[0], ktraj_adc[1], 'r.')  # Plot  sampling points
#plt.show()

# Prepare the sequence output for the scanner
seq.set_definition('FOV', [fovx, fovy, slice_thickness])
seq.set_definition('Name', 'TSE')
seq.set_definition('Nx,Ny', [Nx,Ny])
seq.set_definition('flip', [rf_flip])
seq.set_definition('TR,TE,TI', [TR,TE,TI])
seq.set_definition('SliceThickness', slice_thickness)
seq.set_definition('apodization', apodization)

TEs=int(TE*1000)
TRs=int(TR*1000)
TIs=int(TI*1000)
seq.write(f'IRSE_{rf_flip}TE{TEs:03d}TR{TRs:04d}TI{TIs:03d}.seq')

# Very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within
# slew-rate limits
#rep = seq.test_report()
#print(rep)