# Neuronal Oscillations with Non-sinusoidal Morphology Produce Spurious Phase-to-Amplitude Coupling and Directionality

Code to reproduce simulations published in [Lozano-Soldevilla et al (2016) *Front Comput Neurosci*](https://www.frontiersin.org/articles/10.3389/fncom.2016.00087/full)

Here you will find the scripts to reproduce the figures of the paper (simulations only):

`Figure3A.m`      
`Figure3B.m`          
`Figure4A.m`      
`Figure4B.m`       
`Figure4C.m`      
`Figure4D.m`      
`Figure5.m`      
`bluewhitered.m`

and the core functions to compute bicoherence, cross-frequency coherence and cross-frequency directionality:

`bicoh.m`        
`cfcoh.m`      
`cfd.m`      

Before running the code, you need to first download [fieldtrip](https://github.com/fieldtrip/fieldtrip.git) and add it to your path as follows:

  `restoredefaultpath`        
    `addpath path_to_directory/fieldtrip-master`        
    `ft_defaults`        

In case you have any questions, please do not hesitate to contact me (diegols [arroba] protonmail [punto] ch)
