#BSUB -q 7702ib 
#BSUB -J pro
#BSUB -n 16
#BSUB -o outwater
#BSUB -e errwater
#export OMP_NUM_THREADS=1
#source /fs00/software/lsf/misc/ompthreads.sh

module load gcc/8.3.0
module load fftw/3.3.7-iccifort-17.0.6-avx2
module load zlib/1.2.11
module load impi/5.0.3.048

export PLUMED_KERNEL=/fsa/home/ww_guanxy/enableplumed/plumed/lib/libplumedKernel.so
source /fsa/home/ww_guanxy/enableplumed/gromacs/bin/GMXRC.bash
source /fsa/home/ww_guanxy/amber20/amber.sh
export LD_LIBRARY_PATH=/fsa/home/ww_guanxy/enableplumed/plumed/lib:$LD_LIBRARY_PATH
PATH=/fsa/home/ww_guanxy/amber20/bin:/fsa/home/ww_guanxy/enableplumed/gromacs/bin:/fsa/home/ww_guanxy/anaconda3/envs/py27/bin:/fsa/home/ww_guanxy/amber20/bin:/fsa/home/ww_guanxy/anaconda3/bin:/fs00/software/intel/ps2015u3/impi/5.0.3.048/intel64/bin:/fs00/software/fftw/3.3.7-iccifort-17.0.6-avx2/bin:/fs00/software/gcc/8.3.0/bin:/fsa/home/ww_guanxy/anaconda3/bin:/fsa/home/ww_guanxy/anaconda3/condabin:/fs00/software/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/fs00/software/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/fs00/software/modules/5.0.1/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/fs00/software/bin:/usr/lpp/mmfs/bin:/fsa/home/ww_guanxy/.local/bin:/fsa/home/ww_guanxy/bin

conda activate py27


bash ./tiny0.sh &
bash ./tiny1.sh &
bash ./tiny2.sh &
bash ./tiny3.sh &
bash ./tiny4.sh &
bash ./tiny5.sh &
bash ./tiny6.sh &
bash ./tiny7.sh &
bash ./tiny8.sh &
bash ./tiny9.sh &
bash ./tiny10.sh &
bash ./tiny11.sh &
bash ./tiny12.sh &
bash ./tiny13.sh &
bash ./tiny14.sh &
bash ./tiny15.sh 
