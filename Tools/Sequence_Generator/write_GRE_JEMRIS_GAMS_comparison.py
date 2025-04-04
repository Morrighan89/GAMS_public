"""
GRE sequence (template used for JEMRIS-GAMS comparison in Sub-section 3.2 and for the generation of T1-weighted image in Sub-section 4.2 of https://dx.doi.org/10.2139/ssrn.4882851)

---
"""

from math import pi

import numpy as np
from matplotlib import pyplot as plt
from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

"""
USER INPUTS

These parameters are typically on the user interface of the scanner computer console 
"""
nsa = 1  # Number of averages
n_slices = 1  # Number of slices
Nx = 256
Ny = 256
fovx = 220.e-3  # mm
fovy = 220.e-3
slice_thickness = 2.e-3  # s
slice_gap = 0.e-3  # s
rf_flip = 60  # degrees
rf_offset = 0
apodization=2
print('User inputs setup')


""" SYSTEM LIMITS
Set the hardware limits and initialize sequence object
"""
system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', 
              rf_ringdown_time=10e-6, grad_raster_time=10.e-6,
              rf_dead_time=100e-6)
seq = Sequence(system)

""" TIME CONSTANTS """

TE = np.array([8.e-3])
TR = 500e-3

""" RF

Provide flip angle in deg, first line converts to radians

The function gives the RF sync, the slice selection gradient gz and the rephasing gradient gzr
"""

flip_ang = round(rf_flip * pi / 180, 3)

rf, gz, gzr = make_sinc_pulse(flip_angle=flip_ang, system=system, duration=1.5e-3, slice_thickness=slice_thickness, apodization=apodization,  return_gz=True, time_bw_product=4)
""" READOUT 
Readout gradients and related events
"""

delta_kx = 1 / fovx
gx = make_trapezoid(channel='x', flat_area=Nx * delta_kx, flat_time=5.12e-3, system=system)
adc = make_adc(num_samples=Nx, duration=gx.flat_time, delay=gx.rise_time, system=system)

""" PREPHASE AND REPHASE
"""
delta_ky = 1 / fovy
phase_areas = (np.arange(Ny) - Ny / 2) * delta_ky

gx_pre = make_trapezoid(channel='x', area=-gx.area / 2, duration=3.5e-3, system=system)
gz_reph = make_trapezoid(channel='z', area=-gz.area / 2, duration=3.5e-3, system=system)

""" DELAYS
Echo time (TE) and repetition time (TR). Here, TE is broken down into `delay1` and `delay2`.
"""

delay_TE = TE - calc_duration(gx_pre) - gz.fall_time - gz.flat_time / 2 - calc_duration(
    gx) / 2
delay_TR = TR - calc_duration(gz) - calc_duration(gx_pre) - calc_duration(
    gx) - delay_TE

assert np.all(delay_TE >= 0)
#assert np.all(delay_TR >= calc_duration(gx_spoil, gz_spoil))

""" CONSTRUCT SEQUENCE
"""

rf_phase = 0
rf_inc = 0
freq_offset = gz.amplitude * slice_gap 
rf.freq_offset = freq_offset
seq.add_block(rf, gz)
seq.add_block(make_delay(delay_TR[0]))
for i in range(Ny):
    for j in range(len(TE)):
        rf.phase_offset = rf_phase / 180 * np.pi
        adc.phase_offset = rf_phase / 180 * np.pi

        seq.add_block(rf, gz)
        gy_pre = make_trapezoid(channel='y', area=phase_areas[i], duration=calc_duration(gx_pre), system=system)
        seq.add_block(gx_pre, gy_pre, gz_reph)
        seq.add_block(make_delay(delay_TE[j]))
        seq.add_block(gx, adc)
        gy_pre.amplitude = -gy_pre.amplitude
        seq.add_block(make_delay(delay_TR[j]))
ok, error_report = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print('Timing check passed successfully')
else:
    print('Timing check failed. Error listing follows:')
    [print(e) for e in error_report]

""" PLOTTING TIMNG DIAGRAM """

seq.plot(time_range=(delay_TR[0], delay_TR[0]+0.025))

""" GENERATING `.SEQ` FILE """
# Prepare the sequence output for the scanner
seq.set_definition('FOV', [fovx, fovy, slice_thickness])
seq.set_definition('Name', 'GRE')
seq.set_definition('Nx,Ny', [Nx,Ny])
seq.set_definition('flip', [rf_flip])
seq.set_definition('TR,TE', [TR,TE])
seq.set_definition('apodization', apodization)


seq.write(f'GRE_{rf_flip}TE{int(TE*1000):03d}TR{int(TR*1000):04d}_test.seq')
