"""



Spin Echo (SE) acquisition using the `pypulseq` library. This sequence is typically used for T<sub>2</sub> weighted imaging.
A 2D Fourier transform can be used to reconstruct images from this acquisition. Read more about SE [here](http://mriquestions.com/se-vs-multi-se-vs-fse.html).

---

"""


from math import pi

import numpy as np
import math
import os
from tqdm import tqdm
from matplotlib import pyplot as plt
from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts



def write_SE_basic_sequence(TE,TR,NX,NY,dz,flip):
    """ USER INPUTS
    These parameters are typically on the user interface of the scanner computer console 
    """
    xyanlge=0
    xy_ang = round(xyanlge * pi / 180, 3)
    nsa = 1  # Number of averages
    n_slices = 1  # Number of slices
    Nx = NX
    Ny = NY
    fovx = 304.e-3  # mm
    fovy = 152.e-3
    slice_thickness =dz  # s
    slice_gap = 0.e-3  # s
    rf_flip = flip # degrees
    rf_offset = 0
    apodization=1
    #print('User inputs setup')

    """ SYSTEM LIMITS
    Set the hardware limits and initialize sequence object
    """

    system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', 
                  rf_ringdown_time=10e-6, grad_raster_time=10.e-6,
                  rf_dead_time=100e-6)
    seq = Sequence(system)

    """ TIME CONSTANTS """

    #TE = 500.e-3  # s These are now received as parameter inputs of the function
    #TR = 5  # s
    tau = TE / 2  # s
    readout_time = 4.5e-3
    pre_time = 8e-4  # s

    """ RF """

    flip90 = round(rf_flip * pi / 180, 3) #round to three decimal digits
    flip180 = 180 * pi / 180
    rf90, gz90, _ = make_sinc_pulse(flip_angle=flip90, system=system, duration=2.51e-3, ##Calculate amplitude from flip=360 gamma \int_{t}^{t+duration}B1(t) dt
                                    slice_thickness=slice_thickness, apodization=apodization, 
                                    return_gz=True,time_bw_product=2)
    rf180, gz180, _ = make_sinc_pulse(flip_angle=flip180, system=system, return_gz=True,
                                      duration=2.51e-3, 
                                      slice_thickness=slice_thickness, 
                                      apodization=apodization, 
                                    time_bw_product=2)

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
    #print(gx)

    """ PREPHASE AND REPHASE """

    phase_areas = (np.arange(Ny) - (Ny / 2)) * delta_ky
    gz_reph = make_trapezoid(channel='z', system=system, area=-gz90.area / 2,
                             duration=readout_time / 2)
    gx_pre = make_trapezoid(channel='x', system=system, flat_area=k_width_x / 2, 
                            flat_time=readout_time / 2)
    gy_pre = make_trapezoid(channel='y', system=system, area=phase_areas[-1], 
                            duration=readout_time / 2)


    """ SPOILER """

    #gz_spoil = make_trapezoid(channel='z', system=system, area=gz90.area * 4,
    #                          duration=pre_time * 8)

    """
    Echo time (TE) and repetition time (TR). Here, TE is broken down into `delay1` and `delay2`.
    """

    delay1 = tau - calc_duration(rf90) / 2 - calc_duration(gx_pre)
    #delay1 -= calc_duration(gz_spoil) - calc_duration(rf180) / 2
    delay1 -= calc_duration(rf180) / 2
    delay1 = make_delay(delay1)
    #delay2 = tau - calc_duration(rf180) / 2 - calc_duration(gz_spoil)
    delay2 = tau - calc_duration(rf180) / 2
    delay2 -= calc_duration(gx) / 2
    delay2 = make_delay(delay2)
    delay_TR = TR - calc_duration(rf90) / 2 - calc_duration(gx) / 2 - TE
    #delay_TR -= calc_duration(gy_pre)
    delay_TR = make_delay(delay_TR)
    #print(f'delay_1: {delay1}')
    #print(f'delay_2: {delay2}')
    #print(f'delay_TR: {delay_TR}')

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
            seq.add_block(rf90, gz90)
            gy_pre = make_trapezoid(channel='y', system=system, 
                                    area=phase_areas[i], duration=readout_time / 2)
            seq.add_block(gx_pre, gy_pre, gz_reph)
            seq.add_block(delay1)
            seq.add_block(rf180, gz180)
            seq.add_block(delay2)
            seq.add_block(gx, adc)

            seq.add_block(delay_TR)

    """ PLOTTING TIMING DIAGRAM """

    #seq.plot(time_range=(0.0,0.05))

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
    seq.set_definition('Name', 'SE')
    seq.set_definition('Nx,Ny', [Nx,Ny])
    seq.set_definition('flip', [rf_flip])
    seq.set_definition('TR,TE', [TR,TE])
    seq.set_definition('apodization', apodization)

    filename_seq=f'SE_{rf_flip}TE{int(TE*1000):03d}TR{int(TR*1000):04d}_ap{apodization}_dz{slice_thickness*1000:.1f}_N{Ny:d}'.replace('.', '-') ## USe a similar naming convention for the t2 fit automation tool
    seq.write(os.path.join('test',filename_seq))

    # Very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within
    # slew-rate limits
    #rep = seq.test_report()
    #print(rep)

def main():
    TEs=[10.e-3,20.e-3,50.e-3,120.e-3,320.e-3,500.e-3]
    TRs=[5.0,2.5]
    NXs=[256,128]
    NYs=[128,64]
    dzs=[1.05e-3,2.5e-3,5.e-3,7.5e-3,10.e-3]
    flips=[87]
 
    total_iterations=len(TEs)*len(TRs)*len(NXs)*len(dzs)*len(flips)
    with tqdm(total=total_iterations, desc="Overall Progress") as pbar:
      for TE in TEs:
          for TR in TRs:
             for nx,ny in zip(NXs,NYs):
                for dz in dzs:
                   for flip in flips:
                    write_SE_basic_sequence(TE,TR,nx,ny,dz,flip)
                    pbar.update(1)
    

   


if __name__ == "__main__":
    main()