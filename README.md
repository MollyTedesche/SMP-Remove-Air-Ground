# SMP-Remove-Air-Ground
This is a Matlab Script Module for Systematically Removing Top Air Space and Bottom Ground Space from SnowMicroPenetrometer (SMP) Raw Force Profiles using statistical approach based on the Satyawali et al. (2009) paper: 
"Preliminary characterization of Alpine snow using SnowMicroPen"

PART (1) - Algorithm to remove top of SMP Profile that indicates instrument moving through air space before reaching the snow surface (top 10 mm) 

Rule 1 - first 10 mm of the SMP profile is marked as air because cone tip of the SMP travels about this length in air before hitting the snow surface. Then estimate MPF, SD, & CV values averaged over 10 mm profile length (Mean Force, Standard Deviation, Coefficient of Variation)

Rule 2 - from 11 mm onward till the end of SMP profile, if statistics (MPF,SD,CV) of msmnt are within the range of ±10% of the first 10 mm, it's marked as air. The first set of 10, 1mm consecutive non-air measurements is set as snow surface (snow-air interface is 10 mm thick)

PART (2) - Algorithm to remove bottom of SMP Profile that indicates below snow-ground interface

Rule 1 - Mark the last/bottom 10 mm of the SMP profile as ground, then estimate MPF, SD, & CV values averaged over last 10 mm profile length (Mean Force, Standard Dev, Coeff of Variation)

Rule 2 - From 11 mm from the bottom upward, if statistics (MPF,SD,CV) over each ascending 1 mm length are within the range of ±10% of the last 10mm, it's marked as ground. The first set of 10, 1 mm ascending consecutive non-ground measurements (starting from the bottom up) is set as bottom of the snowpack (snow-ground interface as 10mm thick)

Note that in Satyawali et al., 2009 - the densest/hardest snow types of "Melt-Freeze" and "Hard Crust" have averages of MPF = 6.5N, SD = 1.1N, and CV = 0.24; and maximums of MPF = 23.5N, SD = 5.7N, and CV = 2.0  
   
