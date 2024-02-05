KPL/MK

Description:
   Meta-kernel for n-body orbital propagator

 Author:
   Name: ALESSANDRO 
   Surname: MORSELLI
   Research group: DART
   Department: DAER
   University: Politecnico di Milano 
   Creation: 16/10/2023
   Contact: alessandro.morselli@polimi.it
   Copyright: (c) 2023 A. Morselli, Politecnico di Milano. 
                  All rights reserved.

 Notes:
   This material was prepared to support the course 'Satellite Guidance
   and Navigation', AY 2023/2024.


 NB: this kernel was generated on Windows PC, file paths and line-ending 
     shall be changed on MacOS and linux.

\begindata

    PATH_VALUES = ( 'kernels' )
    PATH_SYMBOLS = ( 'KERNELS' )
    KERNELS_TO_LOAD = (
                       '$KERNELS\naif0012.tls',
                       '$KERNELS\de432s.bsp',
                       '$KERNELS\gm_de432.tpc',
                       '$KERNELS\pck00010.tpc',
                       '$KERNELS\20099942_Apophis.bsp'
                      )

\begintext

KERNELS: 
(downloaded from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/)