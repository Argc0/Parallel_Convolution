# MPI - MPI&OpenMP (Parallel Image Convolution)   

A parallel program that applies convolution filters to images and has been optimized with two different parallel programming techniques: MPI and MPI&OpenMP. This program was implemented by **Ritsogianni Argyro** and **Leonidas Triantafyllou** for the course of *Parallel Systems* (Fall 2017-2018).    

**Example how to run (using terminal):**  

**_MPI:_**  
*mpiexec -n 4 mpi_exe -f waterfall_grey_1920_2520.raw -r 2520 -c 1920 -t grey -s 1 -i 1*  
**_MPI+OPENMP:_**  
*mpiexec -n 4 openmp_exe -f waterfall_1920_2520.raw -r 2520 -c 1920 -t rgb -s 1 -i 1 -d 1*  

# Examples of photos with a slight blurring:    
**_Grey:_**  
![alt text](https://github.com/Argc0/Parallel_Convolution/blob/master/Readme%20and%20photos/grey_images.png)    
**_Rgb:_**  
![alt text](https://github.com/Argc0/Parallel_Convolution/blob/master/Readme%20and%20photos/rgb_images.png)  

# Team Members and Contact Details:  
  * Leonidas Triantafyllou: sdi1400202@di.uoa.gr
  * Argyro Ritsogianni: sdi1400171@di.uoa.gr
