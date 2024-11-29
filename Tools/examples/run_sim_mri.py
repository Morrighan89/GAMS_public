import subprocess
from tqdm import tqdm

def main():
   TEs=[10.e-3,20.e-3,50.e-3,120.e-3,320.e-3,500.e-3]
   TRs=[5.0,2.5]
   NXs=[256,128]
   NYs=[128,64]
   dzs=[1.05e-3,2.5e-3,5.e-3,7.5e-3,10.e-3]
   flips=[90,85]

   total_iterations=len(TEs)*len(TRs)*len(NXs)*len(dzs)*len(flips)
   with tqdm(total=total_iterations, desc="Overall Progress") as pbar:
      for TE in TEs:
          for TR in TRs:
             for nx,ny in zip(NXs,NYs):
                for dz in dzs:
                   for flip in flips:
                    command=f'../../bloch20-11_11-1_ICE4_dp.out NewTub SE_{flip}TE{int(TE*1000):03d}TR{int(TR*1000):04d}_ap1_dz{dz*1000:.1f}_{ny} ../../NewTub4.msh'
                    process=subprocess.run(command)
                    pbar.update(1)


if __name__ == "__main__":
    main()
