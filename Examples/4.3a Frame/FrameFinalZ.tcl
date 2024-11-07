#
#   ____ _ ____         ___                   ____                   ___       _             __                
#  / ___(_)  _ \   _   / _ \ _ __   ___ _ __ / ___|  ___  ___  ___  |_ _|_ __ | |_ ___ _ __ / _| __ _  ___ ___ 
# | |  _| | | | |_| |_| | | | '_ \ / _ \ '_ \\___ \ / _ \/ _ \/ __|  | || '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
# | |_| | | |_| |_   _| |_| | |_) |  __/ | | |___) |  __/  __/\__ \  | || | | | ||  __/ |  |  _| (_| | (_|  __/
#  \____|_|____/  |_|  \___/| .__/ \___|_| |_|____/ \___|\___||___/ |___|_| |_|\__\___|_|  |_|  \__,_|\___\___|
#                           |_|                                                                                
#
# GiD + OpenSees Interface - An Integrated FEA Platform
# Copyright (C) 2016-2020
#
# Lab of R/C and Masonry Structures
# School of Civil Engineering, AUTh
#
# Development Team
#
# T. Kartalis-Kaounis, Dipl. Eng. AUTh, MSc
# V.K. Papanikolaou, Dipl. Eng., MSc DIC, PhD, Asst. Prof. AUTh
#
# Project Contributors
#
# F. Derveni, Dipl. Eng. AUTh
# V.K. Protopapadakis, Dipl. Eng. AUTh, MSc
# T. Papadopoulos, Dipl. Eng. AUTh, MSc
# T. Zachariadis, Dipl. Eng. AUTh, MSc
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

logFile "FrameFinalZ.log"
# --------------------------------------------------------------------------------------------------------------
# U N I T S
# --------------------------------------------------------------------------------------------------------------

# Length : m
# Force  : kN
# Moment : kNm
# Stress : kPa
# Mass   : ton

# --------------------------------------------------------------------------------------------------------------
# M A T E R I A L S / S E C T I O N S ( E L E M E N T S )
# --------------------------------------------------------------------------------------------------------------

# Thermal_Displacement-Based_Beam-Column-C (68)
# Thermal_Displacement-Based_Beam-Column-B (45)

# --------------------------------------------------------------------------------------------------------------
#
# M O D E L  D O M A I N  6DOF  (6)
#
# --------------------------------------------------------------------------------------------------------------

model BasicBuilder -ndm 3 -ndf 6

# --------------------------------------------------------------------------------------------------------------
# N O D E S
# --------------------------------------------------------------------------------------------------------------

# node $NodeTag $XCoord $Ycoord $Zcoord

node      1            0         2.25            0
node      2            0      2.18333            0
node      3            0      2.11667            0
node      4            0         2.05            0
node      5            0      1.98333            0
node      6            0      1.91667            0
node      7            0         1.85            0
node      8            0      1.78333            0
node      9            0      1.71667            0
node     10            0         1.65            0
node     11            0      1.58333            0
node     12            0      1.51667            0
node     13            0         1.45            0
node     14            0      1.38333            0
node     15            0      1.31667            0
node     16            0         1.25            0
node     17            0       1.1875            0
node     18            0        1.125            0
node     19            0        1.125    0.0616667
node     20            0        1.125     0.123333
node     21            0        1.125        0.185
node     22            0        1.125     0.246667
node     23            0        1.125     0.308333
node     24            0        1.125         0.37
node     25            0       1.0625            0
node     26            0        1.125     0.431667
node     27            0        1.125     0.493333
node     28            0            1            0
node     29            0        1.125        0.555
node     30            0        1.125     0.616667
node     31            0        1.125     0.678333
node     32            0     0.933333            0
node     33            0        1.125         0.74
node     34            0        1.125     0.801667
node     35            0     0.866667            0
node     36            0        1.125     0.863333
node     37            0          0.8            0
node     38            0        1.125        0.925
node     39            0        1.125     0.978333
node     40            0     0.733333            0
node     41            0        1.125      1.03167
node     42            0        1.125        1.085
node     43            0     0.666667            0
node     44            0        1.125      1.13833
node     45            0        1.125      1.19167
node     46            0          0.6            0
node     47            0        1.125        1.245
node     48            0     0.533333            0
node     49            0        1.125      1.29833
node     50            0        1.125      1.35167
node     51            0     0.466667            0
node     52            0        1.125        1.405
node     53            0        1.125      1.45833
node     54            0          0.4            0
node     55            0        1.125      1.51167
node     56            0     0.333333            0
node     57            0        1.125        1.565
node     58            0        1.125      1.61833
node     59            0     0.266667            0
node     60            0        1.125      1.67167
node     61            0          0.2            0
node     62            0        1.125        1.725
node     63            0        1.125      1.78667
node     64            0     0.133333            0
node     65            0        1.125      1.84833
node     66            0    0.0666667            0
node     67            0        1.125         1.91
node     68            0            0            0
node     69            0        1.125      1.97167
node     70            0        1.125      2.03333
node     71            0        1.125        2.095
node     72            0        1.125      2.15667
node     73            0        1.125      2.21833
node     74            0        1.125         2.28
node     75            0        1.125      2.34167
node     76            0         2.25         2.65
node     77            0      2.18333         2.65
node     78            0      2.11667         2.65
node     79            0        1.125      2.40333
node     80            0         2.05         2.65
node     81            0      1.98333         2.65
node     82            0      1.91667         2.65
node     83            0         1.85         2.65
node     84            0      1.78333         2.65
node     85            0      1.71667         2.65
node     86            0        1.125        2.465
node     87            0         1.65         2.65
node     88            0      1.58333         2.65
node     89            0      1.51667         2.65
node     90            0        1.125      2.52667
node     91            0         1.45         2.65
node     92            0      1.38333         2.65
node     93            0      1.31667         2.65
node     94            0        1.125      2.58833
node     95            0         1.25         2.65
node     96            0       1.1875         2.65
node     97            0        1.125         2.65
node     98            0       1.0625         2.65
node     99            0            1         2.65
node    100            0     0.933333         2.65
node    101            0     0.866667         2.65
node    102            0          0.8         2.65
node    103            0     0.733333         2.65
node    104            0     0.666667         2.65
node    105            0          0.6         2.65
node    106            0     0.533333         2.65
node    107            0     0.466667         2.65
node    108            0          0.4         2.65
node    109            0     0.333333         2.65
node    110            0     0.266667         2.65
node    111            0          0.2         2.65
node    112            0     0.133333         2.65
node    113            0    0.0666667         2.65
node    114            0            0         2.65

# mass      1            0.01  0.01  0.01  0  0  0
# mass      2            0.01  0.01  0.01  0  0  0
# mass      3            0.01  0.01  0.01  0  0  0
# mass      4            0.01  0.01  0.01  0  0  0
# mass      5            0.01  0.01  0.01  0  0  0
# mass      6            0.01  0.01  0.01  0  0  0
# mass      7            0.01  0.01  0.01  0  0  0
# mass      8            0.01  0.01  0.01  0  0  0
# mass      9            0.01  0.01  0.01  0  0  0
# mass     10            0.01  0.01  0.01  0  0  0
# mass     11            0.01  0.01  0.01  0  0  0
# mass     12            0.01  0.01  0.01  0  0  0
# mass     13            0.01  0.01  0.01  0  0  0
# mass     14            0.01  0.01  0.01  0  0  0
# mass     15            0.01  0.01  0.01  0  0  0
# mass     16            0.01  0.01  0.01  0  0  0
# mass     17            0.01  0.01  0.01  0  0  0
# mass     18            0.01  0.01  0.01  0  0  0
# mass     19            0.01  0.01  0.01  0  0  0
# mass     20            0.01  0.01  0.01  0  0  0
# mass     21            0.01  0.01  0.01  0  0  0
# mass     22            0.01  0.01  0.01  0  0  0
# mass     23            0.01  0.01  0.01  0  0  0
# mass     24            0.01  0.01  0.01  0  0  0
# mass     25            0.01  0.01  0.01  0  0  0
# mass     26            0.01  0.01  0.01  0  0  0
# mass     27            0.01  0.01  0.01  0  0  0
# mass     28            0.01  0.01  0.01  0  0  0
# mass     29            0.01  0.01  0.01  0  0  0
# mass     30            0.01  0.01  0.01  0  0  0
# mass     31            0.01  0.01  0.01  0  0  0
# mass     32            0.01  0.01  0.01  0  0  0
# mass     33            0.01  0.01  0.01  0  0  0
# mass     34            0.01  0.01  0.01  0  0  0
# mass     35            0.01  0.01  0.01  0  0  0
# mass     36            0.01  0.01  0.01  0  0  0
# mass     37            0.01  0.01  0.01  0  0  0
# mass     38            0.01  0.01  0.01  0  0  0
# mass     39            0.01  0.01  0.01  0  0  0
# mass     40            0.01  0.01  0.01  0  0  0
# mass     41            0.01  0.01  0.01  0  0  0
# mass     42            0.01  0.01  0.01  0  0  0
# mass     43            0.01  0.01  0.01  0  0  0
# mass     44            0.01  0.01  0.01  0  0  0
# mass     45            0.01  0.01  0.01  0  0  0
# mass     46            0.01  0.01  0.01  0  0  0
# mass     47            0.01  0.01  0.01  0  0  0
# mass     48            0.01  0.01  0.01  0  0  0
# mass     49            0.01  0.01  0.01  0  0  0
# mass     50            0.01  0.01  0.01  0  0  0
# mass     51            0.01  0.01  0.01  0  0  0
# mass     52            0.01  0.01  0.01  0  0  0
# mass     53            0.01  0.01  0.01  0  0  0
# mass     54            0.01  0.01  0.01  0  0  0
# mass     55            0.01  0.01  0.01  0  0  0
# mass     56            0.01  0.01  0.01  0  0  0
# mass     57            0.01  0.01  0.01  0  0  0
# mass     58            0.01  0.01  0.01  0  0  0
# mass     59            0.01  0.01  0.01  0  0  0
# mass     60            0.01  0.01  0.01  0  0  0
# mass     61            0.01  0.01  0.01  0  0  0
# mass     62            0.01  0.01  0.01  0  0  0
# mass     63            0.01  0.01  0.01  0  0  0
# mass     64            0.01  0.01  0.01  0  0  0
# mass     65            0.01  0.01  0.01  0  0  0
# mass     66            0.01  0.01  0.01  0  0  0
# mass     67            0.01  0.01  0.01  0  0  0
# mass     68            0.01  0.01  0.01  0  0  0
# mass     69            0.01  0.01  0.01  0  0  0
# mass     70            0.01  0.01  0.01  0  0  0
# mass     71            0.01  0.01  0.01  0  0  0
# mass     72            0.01  0.01  0.01  0  0  0
# mass     73            0.01  0.01  0.01  0  0  0
# mass     74            0.01  0.01  0.01  0  0  0
# mass     75            0.01  0.01  0.01  0  0  0
# mass     76            0.01  0.01  0.01  0  0  0
# mass     77            0.01  0.01  0.01  0  0  0
# mass     78            0.01  0.01  0.01  0  0  0
# mass     79            0.01  0.01  0.01  0  0  0
# mass     80            0.01  0.01  0.01  0  0  0
# mass     81            0.01  0.01  0.01  0  0  0
# mass     82            0.01  0.01  0.01  0  0  0
# mass     83            0.01  0.01  0.01  0  0  0
# mass     84            0.01  0.01  0.01  0  0  0
# mass     85            0.01  0.01  0.01  0  0  0
# mass     86            0.01  0.01  0.01  0  0  0
# mass     87            0.01  0.01  0.01  0  0  0
# mass     88            0.01  0.01  0.01  0  0  0
# mass     89            0.01  0.01  0.01  0  0  0
# mass     90            0.01  0.01  0.01  0  0  0
# mass     91            0.01  0.01  0.01  0  0  0
# mass     92            0.01  0.01  0.01  0  0  0
# mass     93            0.01  0.01  0.01  0  0  0
# mass     94            0.01  0.01  0.01  0  0  0
# mass     95            0.01  0.01  0.01  0  0  0
# mass     96            0.01  0.01  0.01  0  0  0
# mass     97            0.01  0.01  0.01  0  0  0
# mass     98            0.01  0.01  0.01  0  0  0
# mass     99            0.01  0.01  0.01  0  0  0
# mass    100            0.01  0.01  0.01  0  0  0
# mass    101            0.01  0.01  0.01  0  0  0
# mass    102            0.01  0.01  0.01  0  0  0
# mass    103            0.01  0.01  0.01  0  0  0
# mass    104            0.01  0.01  0.01  0  0  0
# mass    105            0.01  0.01  0.01  0  0  0
# mass    106            0.01  0.01  0.01  0  0  0
# mass    107            0.01  0.01  0.01  0  0  0
# mass    108            0.01  0.01  0.01  0  0  0
# mass    109            0.01  0.01  0.01  0  0  0
# mass    110            0.01  0.01  0.01  0  0  0
# mass    111            0.01  0.01  0.01  0  0  0
# mass    112            0.01  0.01  0.01  0  0  0
# mass    113            0.01  0.01  0.01  0  0  0
# mass    114            0.01  0.01  0.01  0  0  0




# --------------------------------------------------------------------------------------------------------------
# R E S T R A I N T S
# --------------------------------------------------------------------------------------------------------------

# fix $NodeTag x-transl y-transl z-transl x-rot y-rot z-rot

fix      1   1   0   1   0   1   1
fix     68   1   1   1   0   1   1
fix     76   1   0   1   0   1   1
fix    114   1   1   1   0   1   1

fix    59    1   0   1   0   0  0
fix    61    1   0   1   0   0  0
fix    64    1   0   1   0   0  0
fix    66    1   0   1   0   0  0
fix    110   1   0   1   0  0  0
fix    111   1   0   1   0  0  0
fix    112   1   0   1   0  0  0
fix    113   1   0   1   0  0  0
#fix     18   0   1   1   1   1   1
# fix     97   0   0   1   1   0   0
#fix     17   0   0   1   1   0   0
#fix     96   0   0   1   1   0   0

# --------------------------------------------------------------------------------------------------------------
# D I S P L A C E M E N T - B A S E D   B E A M - C O L U M N   E L E M E N T S
# --------------------------------------------------------------------------------------------------------------

# Geometric Transformation
geomTransf Linear 1  -1.0 0.0 0; #Vertical with angle = 0
geomTransf Linear 2  0 0 1; #Non-vertical with angle = 0
geomTransf PDelta 3  -1.0 0.0 0; #Vertical with angle = 0
geomTransf PDelta 4  0 0 1; #Non-vertical with angle = 0
geomTransf Corotational 5  -1.0 0.0 0; #Vertical with angle = 0
geomTransf Corotational 6  0 0 1; #Non-vertical with angle = 0


# Sections Definition used by dispBeamColumn Elements
# (if they have not already been defined on this model domain)

######### Section for element: Thermal_Displacement-Based_Beam-Column-C  ##########
uniaxialMaterial DamagePlasticityConcreteECT 618 -36300 -0.0025 -32000 -0.02 0.2 2900 2.72727e+06
uniaxialMaterial DamagePlasticityConcreteECT 658 -36300 -0.0025 -32000 -0.02 0.2 2900 2.72727e+06
uniaxialMaterial DamagePlasticityConcreteECT 678 -60300 -0.0025 -50000 -0.02 0.2 2900 2.72727e+06
uniaxialMaterial DPMsteelEC 668 351000 1.91e+08 0.01 0 1 0 1


section fiberSecThermal 669 -GJ 1e10  {

# Create the Core fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 658     20     20  -0.125000  -0.125000   0.125000   0.125000


# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 668  2   0.00020110   0.085000   0.085000   0.085000  -0.085000
layer straight 668  2   0.00020110  -0.085000   0.085000  -0.085000  -0.085000

# Create the middle bars

# fiber $yLoc $zLoc $A $matTag
fiber   0.085000 0   0.00020110 668
fiber  -0.085000 0   0.00020110 668

# Create the middle Reinforcing Bars along local y axis

# fiber $yLoc $zLoc $A $matTag
fiber 0   0.085000     0.00020110 668
fiber 0  -0.085000     0.00020110 668
}
########### Joint in column ##############
section fiberSecThermal 659 -GJ 1e10  {

# Create the Core fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 678     20     20  -0.125000  -0.125000   0.125000   0.125000


# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 668  2   0.00020110   0.085000   0.085000   0.085000  -0.085000
layer straight 668  2   0.00020110  -0.085000   0.085000  -0.085000  -0.085000

# Create the middle bars

# fiber $yLoc $zLoc $A $matTag
fiber   0.085000 0   0.00020110 668
fiber  -0.085000 0   0.00020110 668

# Create the middle Reinforcing Bars along local y axis

# fiber $yLoc $zLoc $A $matTag
fiber 0   0.085000     0.00020110 668
fiber 0  -0.085000     0.00020110 668
}
################# Section for element: Thermal_Displacement-Based_Beam-Column-B ############
uniaxialMaterial SteelDPM1 667 408000 1.88e+08 0.01 0 1 0 1

section fiberSecThermal 670 -GJ 1e10  {

# Create the Core fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618      16     24  -0.100000  -0.125000   0.100000   0.125000

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 667   2    0.00028350  -0.060000   0.085000   0.060000   0.085000

# Create the Bottom bars (face on local z negative dir)

layer straight 667   2   0.00028350  -0.060000  -0.085000   0.060000  -0.085000
}

########### Joint in Beam ###########
section fiberSecThermal 671 -GJ 1e10  {

# Create the Core fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 678      16     24  -0.100000  -0.125000   0.100000   0.125000

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 667   2    0.00028350  -0.060000   0.085000   0.060000   0.085000

# Create the Bottom bars (face on local z negative dir)

layer straight 667   2   0.00028350  -0.060000  -0.085000   0.060000  -0.085000
}

# Displacement-Based Beam Column Element Definition

# element dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>

element dispBeamColumnThermal      1     68     66  3    669  6   -mass 0.490625
element dispBeamColumnThermal      2     66     64  3    669  6   -mass 0.490625
element dispBeamColumnThermal      3     64     61  3    669  6   -mass 0.490625
element dispBeamColumnThermal      4     61     59  3    669  6   -mass 0.490625
element dispBeamColumnThermal      5     59     56  3    669  6   -mass 0.490625
element dispBeamColumnThermal      6     56     54  3    669  6   -mass 0.490625
element dispBeamColumnThermal      7     54     51  3    669  6   -mass 0.490625
element dispBeamColumnThermal      8     51     48  3    669  6   -mass 0.490625
element dispBeamColumnThermal      9     48     46  3    669  6   -mass 0.490625
element dispBeamColumnThermal     10     46     43  3    669  6   -mass 0.490625
element dispBeamColumnThermal     11     43     40  3    669  6   -mass 0.490625
element dispBeamColumnThermal     12     40     37  3    669  6   -mass 0.490625
element dispBeamColumnThermal     13     37     35  3    669  6   -mass 0.490625
element dispBeamColumnThermal     14     35     32  3    669  6   -mass 0.490625
element dispBeamColumnThermal     15     32     28  3    669  6   -mass 0.490625
element dispBeamColumnThermal     16     28     25  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     17     25     18  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     18     18     17  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     19     17     16  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     20     16     15  3    669  6   -mass 0.490625
element dispBeamColumnThermal     21     15     14  3    669  6   -mass 0.490625
element dispBeamColumnThermal     22     14     13  3    669  6   -mass 0.490625
element dispBeamColumnThermal     23     13     12  3    669  6   -mass 0.490625
element dispBeamColumnThermal     24     12     11  3    669  6   -mass 0.490625
element dispBeamColumnThermal     25     11     10  3    669  6   -mass 0.490625
element dispBeamColumnThermal     26     10      9  3    669  6   -mass 0.490625
element dispBeamColumnThermal     27      9      8  3    669  6   -mass 0.490625
element dispBeamColumnThermal     28      8      7  3    669  6   -mass 0.490625
element dispBeamColumnThermal     29      7      6  3    669  6   -mass 0.490625
element dispBeamColumnThermal     30      6      5  3    669  6   -mass 0.490625
element dispBeamColumnThermal     31      5      4  3    669  6   -mass 0.490625
element dispBeamColumnThermal     32      4      3  3    669  6   -mass 0.490625
element dispBeamColumnThermal     33      3      2  3    669  6   -mass 0.490625
element dispBeamColumnThermal     34      2      1  3    669  6   -mass 0.490625
element dispBeamColumnThermal     35     97     96  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     36     96     95  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     37     95     93  3    669  6   -mass 0.490625
element dispBeamColumnThermal     38     93     92  3    669  6   -mass 0.490625
element dispBeamColumnThermal     39     92     91  3    669  6   -mass 0.490625
element dispBeamColumnThermal     40     91     89  3    669  6   -mass 0.490625
element dispBeamColumnThermal     41     89     88  3    669  6   -mass 0.490625
element dispBeamColumnThermal     42     88     87  3    669  6   -mass 0.490625
element dispBeamColumnThermal     43     87     85  3    669  6   -mass 0.490625
element dispBeamColumnThermal     44     85     84  3    669  6   -mass 0.490625
element dispBeamColumnThermal     45     84     83  3    669  6   -mass 0.490625
element dispBeamColumnThermal     46     83     82  3    669  6   -mass 0.490625
element dispBeamColumnThermal     47     82     81  3    669  6   -mass 0.490625
element dispBeamColumnThermal     48     81     80  3    669  6   -mass 0.490625
element dispBeamColumnThermal     49     80     78  3    669  6   -mass 0.490625
element dispBeamColumnThermal     50     78     77  3    669  6   -mass 0.490625
element dispBeamColumnThermal     51     77     76  3    669  6   -mass 0.490625
element dispBeamColumnThermal     52     97     98  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     53     98     99  3    659  6   -mass 0.490625 #
element dispBeamColumnThermal     54     99    100  3    669  6   -mass 0.490625
element dispBeamColumnThermal     55    100    101  3    669  6   -mass 0.490625
element dispBeamColumnThermal     56    101    102  3    669  6   -mass 0.490625
element dispBeamColumnThermal     57    102    103  3    669  6   -mass 0.490625
element dispBeamColumnThermal     58    103    104  3    669  6   -mass 0.490625
element dispBeamColumnThermal     59    104    105  3    669  6   -mass 0.490625
element dispBeamColumnThermal     60    105    106  3    669  6   -mass 0.490625
element dispBeamColumnThermal     61    106    107  3    669  6   -mass 0.490625
element dispBeamColumnThermal     62    107    108  3    669  6   -mass 0.490625
element dispBeamColumnThermal     63    108    109  3    669  6   -mass 0.490625
element dispBeamColumnThermal     64    109    110  3    669  6   -mass 0.490625
element dispBeamColumnThermal     65    110    111  3    669  6   -mass 0.490625
element dispBeamColumnThermal     66    111    112  3    669  6   -mass 0.490625
element dispBeamColumnThermal     67    112    113  3    669  6   -mass 0.490625
element dispBeamColumnThermal     68    113    114  3    669  6   -mass 0.490625
element dispBeamColumnThermal     69     18     19  3    671  5   -mass   0.3925 #
element dispBeamColumnThermal     70     19     20  3    670  5   -mass   0.3925
element dispBeamColumnThermal     71     20     21  3    670  5   -mass   0.3925
element dispBeamColumnThermal     72     21     22  3    670  5   -mass   0.3925
element dispBeamColumnThermal     73     22     23  3    670  5   -mass   0.3925
element dispBeamColumnThermal     74     23     24  3    670  5   -mass   0.3925
element dispBeamColumnThermal     75     24     26  3    670  5   -mass   0.3925
element dispBeamColumnThermal     76     26     27  3    670  5   -mass   0.3925
element dispBeamColumnThermal     77     27     29  3    670  5   -mass   0.3925
element dispBeamColumnThermal     78     29     30  3    670  5   -mass   0.3925
element dispBeamColumnThermal     79     30     31  3    670  5   -mass   0.3925
element dispBeamColumnThermal     80     31     33  3    670  5   -mass   0.3925
element dispBeamColumnThermal     81     33     34  3    670  5   -mass   0.3925
element dispBeamColumnThermal     82     34     36  3    670  5   -mass   0.3925
element dispBeamColumnThermal     83     36     38  3    670  5   -mass   0.3925
element dispBeamColumnThermal     84     38     39  3    670  5   -mass   0.3925
element dispBeamColumnThermal     85     39     41  3    670  5   -mass   0.3925
element dispBeamColumnThermal     86     41     42  3    670  5   -mass   0.3925
element dispBeamColumnThermal     87     42     44  3    670  5   -mass   0.3925
element dispBeamColumnThermal     88     44     45  3    670  5   -mass   0.3925
element dispBeamColumnThermal     89     45     47  3    670  5   -mass   0.3925
element dispBeamColumnThermal     90     47     49  3    670  5   -mass   0.3925
element dispBeamColumnThermal     91     49     50  3    670  5   -mass   0.3925
element dispBeamColumnThermal     92     50     52  3    670  5   -mass   0.3925
element dispBeamColumnThermal     93     52     53  3    670  5   -mass   0.3925
element dispBeamColumnThermal     94     53     55  3    670  5   -mass   0.3925
element dispBeamColumnThermal     95     55     57  3    670  5   -mass   0.3925
element dispBeamColumnThermal     96     57     58  3    670  5   -mass   0.3925
element dispBeamColumnThermal     97     58     60  3    670  5   -mass   0.3925
element dispBeamColumnThermal     98     60     62  3    670  5   -mass   0.3925
element dispBeamColumnThermal     99     62     63  3    670  5   -mass   0.3925
element dispBeamColumnThermal    100     63     65  3    670  5   -mass   0.3925
element dispBeamColumnThermal    101     65     67  3    670  5   -mass   0.3925
element dispBeamColumnThermal    102     67     69  3    670  5   -mass   0.3925
element dispBeamColumnThermal    103     69     70  3    670  5   -mass   0.3925
element dispBeamColumnThermal    104     70     71  3    670  5   -mass   0.3925
element dispBeamColumnThermal    105     71     72  3    670  5   -mass   0.3925
element dispBeamColumnThermal    106     72     73  3    670  5   -mass   0.3925
element dispBeamColumnThermal    107     73     74  3    670  5   -mass   0.3925
element dispBeamColumnThermal    108     74     75  3    670  5   -mass   0.3925
element dispBeamColumnThermal    109     75     79  3    670  5   -mass   0.3925
element dispBeamColumnThermal    110     79     86  3    670  5   -mass   0.3925
element dispBeamColumnThermal    111     86     90  3    670  5   -mass   0.3925
element dispBeamColumnThermal    112     90     94  3    670  5   -mass   0.3925
element dispBeamColumnThermal    113     94     97  3    671  5   -mass   0.3925 #

# --------------------------------------------------------------------------------------------------------------
#
# D O M A I N  C O M M O N S
#
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# R E C O R D E R S
# --------------------------------------------------------------------------------------------------------------

recorder Node -file Node_displacements.out -time -nodeRange 1 114 -dof 1 2 3 disp
recorder Node -file Node_rotations.out -time -nodeRange 1 114 -dof 4 5 6 disp
recorder Node -file Node_forceReactions.out -time -nodeRange 1 114 -dof 1 2 3 reaction
recorder Node -file Node_momentReactions.out -time -nodeRange 1 114 -dof 4 5 6 reaction
recorder Element -file DispBeamColumn_localForce.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 localForce
recorder Element -file DispBeamColumn_basicDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 basicDeformation
recorder Element -file DispBeamColumn_plasticDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 plasticDeformation

#recorder Element -file fiberUP.out -time -ele 92 section 1 fiber 0.085000 0.085000 667 stressStrain

puts ""
puts " __   __       __          __                   _       "
puts "/ _ .|  \\ _|_ /  \\ _  _ _ (_  _ _ _  | _ |_ _ _(_ _  _ _"
puts "\\__)||__/  |  \\__/|_)(-| )__)(-(-_)  || )|_(-| | (_|(_(-"
puts "                  |                                     "
puts "                                                  v2.8.0"
puts "Analysis summary"
puts ""
puts "Interval 1 : Static - [expr int(1+100)] steps"
puts "Interval 2 : Static - [expr int(1+360)] steps"
puts ""
puts "----------------"
puts ""
set time_start [clock seconds]
puts "Starting analysis at [clock format $time_start -format %H:%M:%S]\n"

# --------------------------------------------------------------------------------------------------------------
#
# I N T E R V A L   1
#
# --------------------------------------------------------------------------------------------------------------

puts "Running interval 1\n"
# Loads - Plain Pattern
pattern Plain 100 Linear {
    load     38        0    -30.9        0        0        0        0
    load     62        0    -30.9        0        0        0        0
	#load     34        0    -0.5       0        0        0        0
    #load     51        0    -0.5        0        0        0        0
	
	#eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wy $Wz <$Wx>
	
	
	# eleLoad -ele   69  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   70  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   71  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   72  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   73  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   74  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   75  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   76  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   77  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   78  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   79  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   80  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   81  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   82  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   83  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   84  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   85  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   86  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   87  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   88  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   89  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   90  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   91  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   92  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   93  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   94  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   95  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   96  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   97  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   98  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   99  -type -beamUniform   0  -1.3  0
    # eleLoad -ele   100  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   101  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   102  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   103  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   104  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   105  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   106  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   107  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   108  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   109  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   110  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   111  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   112  -type -beamUniform  0  -1.3  0
    # eleLoad -ele   113  -type -beamUniform  0  -1.3  0
    
	# eleLoad -ele    1   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    2   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    3   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    4   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    5   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    6   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    7   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    8   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    9   -type -beamUniform   0  0  -1.6
	# eleLoad -ele    10   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    11   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    12   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    13   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    14   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    15   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    16   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    17   -type -beamUniform  0  0  -1.6
    
	# eleLoad -ele    52  -type -beamUniform   0  0  -1.6
	# eleLoad -ele    53  -type -beamUniform   0  0  -1.6
	# eleLoad -ele    54   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    55   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    56   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    57   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    58   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    59   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    60   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    61   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    62   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    63   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    64   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    65   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    66   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    67   -type -beamUniform  0  0  -1.6
	# eleLoad -ele    68   -type -beamUniform	 0  0  -1.6
}
# recording the initial status

record

# Analysis options

system BandGeneral
numberer RCM
constraints Transformation
integrator LoadControl 0.01
test RelativeEnergyIncr 0.001 50 2
algorithm Newton
analysis Static
set Lincr 0.01
set Nsteps 100
set committedSteps 1
set LoadCounter 0


set strIni {}
variable testTypeStatic RelativeEnergyIncr
variable TolStatic 0.001
variable maxNumIterStatic 50
variable algorithmTypeStatic Newton

for {set i 1} { $i <= $Nsteps } {incr i 1} {
    set t [format "%7.5f" [expr [getTime] + $Lincr]]
    puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
    set AnalOk [analyze 1]
    if {$AnalOk !=0} {
        break
    } else {
        set LoadCounter [expr $LoadCounter+1.0]
        set committedSteps [expr $committedSteps+1]
    }
}

if {$AnalOk != 0} {; # if analysis fails, alternative algorithms and substepping is applied
    set firstFail 1
    set AnalOk 0
    set Nk 1
    set returnToInitStepFlag 0
    while {$LoadCounter < $Nsteps && $AnalOk == 0} {
        if {($Nk==2.0 && $AnalOk==0) || ($Nk==1 && $AnalOk==0)} {
            set Nk 1
            if {$returnToInitStepFlag} {
                puts "\nBack to initial step\n"
                set returnToInitStepFlag 0
            }
            if {$firstFail == 0} { # for the first time only, do not repeat previous failed step
                integrator LoadControl $Lincr; # reset to original increment
                set t [format "%7.5f" [expr [getTime] + $Lincr]]
                puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
                set AnalOk [analyze 1]; # zero for convergence
            } else {
                set AnalOk 1
                set firstFail 0
            }
            if {$AnalOk == 0} {
                set LoadCounter [expr $LoadCounter+1.0/$Nk]
                set committedSteps [expr $committedSteps+1]
            }
        }; # end if Nk=1
        # substepping /2
        if {($AnalOk !=0 && $Nk==1) || ($AnalOk==0 && $Nk==4.0)} {
            set Nk 2.0; # reduce step size
            set continueFlag 1
            puts "\nInitial step is divided by 2\n"
            set LincrReduced [expr $Lincr/$Nk]
            integrator LoadControl $LincrReduced
            for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $LincrReduced]]
                puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
                set AnalOk [analyze 1]; # zero for convergence
                if {$AnalOk == 0} {
                    set LoadCounter [expr $LoadCounter+1.0/$Nk]
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=2.0
        # substepping /4
        if {($AnalOk !=0 && $Nk==2.0) || ($AnalOk==0 && $Nk==8.0)} {
            set Nk 4.0; # reduce step size
            set continueFlag 1
            puts "\nInitial step is divided by 4\n"
            set LincrReduced [expr $Lincr/$Nk]
            integrator LoadControl $LincrReduced
            for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $LincrReduced]]
                puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
                set AnalOk [analyze 1]; # zero for convergence
                if {$AnalOk == 0} {
                    set LoadCounter [expr $LoadCounter+1.0/$Nk]
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=4
        # substepping /8
        if {$AnalOk !=0 && $Nk==4.0 || ($Nk == 16.0 && $AnalOk == 0)} {
            set Nk 8.0; # reduce step size
            set continueFlag 1
            puts "\nInitial step is divided by 8\n"
            set LincrReduced [expr $Lincr/$Nk]
            integrator LoadControl $LincrReduced
            for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $LincrReduced]]
                puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
                set AnalOk [analyze 1]; # zero for convergence
                if {$AnalOk == 0} {
                    set LoadCounter [expr $LoadCounter+1.0/$Nk]
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=8
        # substepping /16
        if {($Nk == 8 && $AnalOk!=0)} {
            set Nk 16.0; # reduce step size
            set continueFlag 1
            puts "\nInitial step is divided by 16\n"
            set LincrReduced [expr $Lincr/$Nk]
            integrator LoadControl $LincrReduced
            for {set ik 1} {$ik <=$Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $LincrReduced]]
                puts -nonewline "(1) $algorithmTypeStatic$strIni LF $t "
                set AnalOk [analyze 1]; # zero for convergence
                if {$AnalOk == 0} {
                    set LoadCounter [expr $LoadCounter+1.0/$Nk]
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=16
    }; # end while loop
}; # end if AnalOk

if {$AnalOk == 0} {
    puts "\nAnalysis completed SUCCESSFULLY"
    puts "Committed steps : $committedSteps\n"
} else {
    puts "\nAnalysis FAILED"
    puts "Committed steps : $committedSteps\n"
}

# all previously defined patterns are constant for so on.

loadConst -time 0.0

# --------------------------------------------------------------------------------------------------------------
#
# I N T E R V A L   2
#
# --------------------------------------------------------------------------------------------------------------

puts "Running interval 2\n"

pattern Plain 200 Linear {
    eleLoad -ele   69  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   70  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   71  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   72  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   73  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   74  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   75  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   76  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   77  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   78  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   79  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   80  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   81  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   82  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   83  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   84  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   85  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   86  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   87  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   88  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   89  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   90  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   91  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   92  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   93  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   94  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   95  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   96  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   97  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   98  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   99  -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   100 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   101 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   102 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   103 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   104 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   105 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   106 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   107 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   108 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   109 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   110 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   111 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   112 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
    eleLoad -ele   113 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation   -0.1    0.1     -0.125      0.125       16       24
	#eleLoad -ele    1   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    #eleLoad -ele    2   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    #eleLoad -ele    3   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    #eleLoad -ele    4   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    5   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    6   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
	eleLoad -ele    7   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    8   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    9   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    10   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    11   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    12   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    13   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
	eleLoad -ele    14   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    15   -type -beamThermal -z   -source  "ColumnRC1.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    16   -type -beamThermal -z   -source  "ColumnRC3.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    17   -type -beamThermal -z   -source  "ColumnRC3.dat"  -genInterpolation     -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    52  -type -beamThermal -z   -source    "ColumnRC3.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    53  -type -beamThermal -z   -source    "ColumnRC3.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    54   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    55   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    56   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    57   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    58   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    59   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    60   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    61   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    62   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    63   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    eleLoad -ele    64   -type -beamThermal -z   -source   "ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
	#eleLoad -ele    65   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    #eleLoad -ele    66   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20  
	#eleLoad -ele    67   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20
    #eleLoad -ele    68   -type -beamThermal -z   -source   "../Records/Thermal_Load/ColumnRC1.dat"  -genInterpolation    -0.125      0.125     -0.125      0.125       20      20 
}
# recording the initial status

record

# Analysis options


system UmfPack
numberer RCM
constraints Transformation
integrator Newmark 0.45 0.8
test RelativeEnergyIncr 0.001 500 2
algorithm KrylovNewton
analysis Transient

# Loads - Uniform Excitation

set iGMdirection 1
set DtAnalysis 60
set TmaxAnalysis 10800

set committedSteps 1
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)]

set strIni {}
variable testTypeDynamic RelativeEnergyIncr
variable TolDynamic 0.001
variable maxNumIterDynamic 100
variable algorithmTypeDynamic Newton

for {set i 1} { $i <= $Nsteps } {incr i 1} {
    set t [format "%7.5f" [expr [getTime] + $DtAnalysis]]
    puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
    set AnalOk [analyze 1 $DtAnalysis]; # perform analysis - returns 0 if analysis was successful
    if {$AnalOk == 0} {
        set committedSteps [expr $committedSteps+1]
    } else {
        break
    }
}

if {$AnalOk != 0} {; # if analysis fails, alternative algorithms and substepping is applied
    set firstFail 1
    set AnalOk 0
    set controlTime [getTime]
    set Nk 1; # dt = dt/Nk
    set returnToInitStepFlag 0
    while {$controlTime < $TmaxAnalysis && $AnalOk == 0} {
        if { ($Nk == 1 && $AnalOk == 0) || ($Nk == 2 && $AnalOk == 0) } {
            set Nk 1
            if {$returnToInitStepFlag} {
                puts "\nBack to initial step\n"
                set returnToInitStepFlag 0
            }
            if {$firstFail == 0} {; # for the first time only, do not repeat previous failed step
                set t [format "%7.5f" [expr [getTime] + $DtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $DtAnalysis]
            } else {
                set AnalOk 1
                set firstFail 0
            }
            if {$AnalOk == 0} {
                set committedSteps [expr $committedSteps+1]
            }
        }; # end if Nk=1
        # substepping /2
        if {($Nk == 1 && $AnalOk!=0) || ($Nk == 4 && $AnalOk==0)} {
            set Nk 2.0
            set continueFlag 1
            puts "\nInitial step is divided by 2\n"
            set currTime1 [getTime]
            set curStep [expr int($currTime1/$DtAnalysis)]
            set remStep1 [expr int(($Nsteps-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $DtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=2
        # substepping /4
        if {($Nk == 2 && $AnalOk!=0) || ($Nk == 8 && $AnalOk == 0)} {
            set Nk 4.0
            set continueFlag 1
            puts "\nInitial step is divided by 4\n"
            set currTime2 [getTime]
            set curStep [expr ($currTime2-$currTime1)/$ReducedDtAnalysis]
            set remStep2 [expr int(($remStep1-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=4
        # substepping /8
        if {($Nk == 4 && $AnalOk!=0) || ($Nk == 16 && $AnalOk == 0)} {
            set Nk 8.0
            set continueFlag 1
            puts "\nInitial step is divided by 8\n"
            set currTime3 [getTime]
            set curStep [expr ($currTime3-$currTime2)/$ReducedDtAnalysis]
            set remStep3 [expr int(($remStep2-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=8
        # substepping /16
        if {($Nk == 8 && $AnalOk!=0)} {
            set Nk 16.0
            set continueFlag 1
            puts "\nInitial step is divided by 16\n"
            set currTime4 [getTime]
            set curStep [expr ($currTime4-$currTime3)/$ReducedDtAnalysis]
            set remStep4 [expr int(($remStep3-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=16
        set controlTime [getTime]
    }
}

if {$AnalOk == 0} {
    puts "\nAnalysis completed SUCCESSFULLY"
    puts "Committed steps : $committedSteps\n"
} else {
    puts "\nAnalysis FAILED"
    puts "Committed steps : $committedSteps\n"
}
# all previously defined patterns are constant for so on.

loadConst -time 0.0

# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
#
# I N T E R V A L   3
#
# --------------------------------------------------------------------------------------------------------------

puts "Running interval 3\n"

pattern Plain 300 Linear {
    load     49        0    -25        0        0        0        0
    #load     62        0    -5        0        0        0        0
	
}
# recording the initial status

record

# Analysis options


system UmfPack
numberer RCM
constraints Transformation
integrator Newmark 0.25 0.5
test RelativeEnergyIncr 0.0001 500 2
algorithm KrylovNewton
analysis Transient

# Loads - Uniform Excitation

set iGMdirection 1
set DtAnalysis 0.01
set TmaxAnalysis 1000

set committedSteps 1
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)]

set strIni {}
variable testTypeDynamic RelativeEnergyIncr
variable TolDynamic 0.0001
variable maxNumIterDynamic 500
variable algorithmTypeDynamic Newton

for {set i 1} { $i <= $Nsteps } {incr i 1} {
    set t [format "%7.5f" [expr [getTime] + $DtAnalysis]]
    puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
    set AnalOk [analyze 1 $DtAnalysis]; # perform analysis - returns 0 if analysis was successful
    if {$AnalOk == 0} {
        set committedSteps [expr $committedSteps+1]
    } else {
        break
    }
}

if {$AnalOk != 0} {; # if analysis fails, alternative algorithms and substepping is applied
    set firstFail 1
    set AnalOk 0
    set controlTime [getTime]
    set Nk 1; # dt = dt/Nk
    set returnToInitStepFlag 0
    while {$controlTime < $TmaxAnalysis && $AnalOk == 0} {
        if { ($Nk == 1 && $AnalOk == 0) || ($Nk == 2 && $AnalOk == 0) } {
            set Nk 1
            if {$returnToInitStepFlag} {
                puts "\nBack to initial step\n"
                set returnToInitStepFlag 0
            }
            if {$firstFail == 0} {; # for the first time only, do not repeat previous failed step
                set t [format "%7.5f" [expr [getTime] + $DtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $DtAnalysis]
            } else {
                set AnalOk 1
                set firstFail 0
            }
            if {$AnalOk == 0} {
                set committedSteps [expr $committedSteps+1]
            }
        }; # end if Nk=1
        # substepping /2
        if {($Nk == 1 && $AnalOk!=0) || ($Nk == 4 && $AnalOk==0)} {
            set Nk 2.0
            set continueFlag 1
            puts "\nInitial step is divided by 2\n"
            set currTime1 [getTime]
            set curStep [expr int($currTime1/$DtAnalysis)]
            set remStep1 [expr int(($Nsteps-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $DtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=2
        # substepping /4
        if {($Nk == 2 && $AnalOk!=0) || ($Nk == 8 && $AnalOk == 0)} {
            set Nk 4.0
            set continueFlag 1
            puts "\nInitial step is divided by 4\n"
            set currTime2 [getTime]
            set curStep [expr ($currTime2-$currTime1)/$ReducedDtAnalysis]
            set remStep2 [expr int(($remStep1-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=4
        # substepping /8
        if {($Nk == 4 && $AnalOk!=0) || ($Nk == 16 && $AnalOk == 0)} {
            set Nk 8.0
            set continueFlag 1
            puts "\nInitial step is divided by 8\n"
            set currTime3 [getTime]
            set curStep [expr ($currTime3-$currTime2)/$ReducedDtAnalysis]
            set remStep3 [expr int(($remStep2-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=8
        # substepping /16
        if {($Nk == 8 && $AnalOk!=0)} {
            set Nk 16.0
            set continueFlag 1
            puts "\nInitial step is divided by 16\n"
            set currTime4 [getTime]
            set curStep [expr ($currTime4-$currTime3)/$ReducedDtAnalysis]
            set remStep4 [expr int(($remStep3-$curStep)*2.0)]
            set ReducedDtAnalysis [expr $ReducedDtAnalysis/2.0]
            for {set ik 1} {$ik <= $Nk} {incr ik 1} {
                if {$continueFlag==0} {
                    break
                }
                set t [format "%7.5f" [expr [getTime] + $ReducedDtAnalysis]]
                puts -nonewline "(1) $algorithmTypeDynamic$strIni Time $t "
                set AnalOk [analyze 1 $ReducedDtAnalysis]
                if {$AnalOk == 0} {
                    set committedSteps [expr $committedSteps+1]
                } else {
                    set continueFlag 0
                }
            }
            if {$AnalOk == 0} {
                set returnToInitStepFlag 1
            }
        }; # end if Nk=16
        set controlTime [getTime]
    }
}

if {$AnalOk == 0} {
    puts "\nAnalysis completed SUCCESSFULLY"
    puts "Committed steps : $committedSteps\n"
} else {
    puts "\nAnalysis FAILED"
    puts "Committed steps : $committedSteps\n"
}

# all previously defined patterns are constant for so on.

loadConst -time 0.0
# #####################
set hour 0.0
set minute 0.0
set second 0.0
set time_end [clock seconds]
set analysisTime [expr $time_end-$time_start]

puts "Analysis finished at [clock format $time_end -format %H:%M:%S]\n"

if {$analysisTime<60} {
    if {$analysisTime==0} {
        puts "Analysis time : less than one second"
    } elseif {$analysisTime==1} {
        puts "Analysis time : 1 second"
    } else {
        puts "Analysis time : $analysisTime seconds"
    }

}  elseif {$analysisTime<3600} {
    set minutes [expr $analysisTime/60]
    set seconds [expr $analysisTime%60]

    if {$minutes==1} {
        puts -nonewline "Analysis time : 1 minute"
    } else {
        puts -nonewline "Analysis time : $minutes minutes"
    }

    if {$seconds==0} {
        puts ""
    } elseif {$seconds==1} {
        puts " and 1 second"
    } else {
        puts " and $seconds seconds"
    }

} else  {
    set hours [expr $analysisTime/3600]
    set minutes [expr ($analysisTime%3600)/60]
    set seconds [expr ($analysisTime%3600)%60]

    if {$hours==1} {
        puts -nonewline "Analysis time : 1 hour"
    } else {
        puts -nonewline "Analysis time : $hours hours"
    }

    if {$minutes==0} {
    } elseif {$minute==1} {
        puts -nonewline ", 1 minute"
    } else {
        puts -nonewline ", $minutes minutes"
    }

    if {$seconds==0} {
        puts ""
    } elseif {$second==1} {
        puts " and 1 second"
    } else {
        puts " and $seconds seconds"
    }
}

# --------------------------------------------------------------------------------------------------------------
#
# M E T A D A T A
#
# --------------------------------------------------------------------------------------------------------------

# Number of nodes
# 114

# Elements 1D
# 113

# Elements 2D
# 0

# Elements 3D
# 0

# DispBeamColumn
# 113
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# F R A M E   L O C A L   A X E S   O R I E N T A T I O N
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#      ID                           Type                       Local-x                       Local-y                       Local-z          Literal      Material / Section
#
#       1                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       2                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       3                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       4                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       5                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       6                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       7                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       8                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#       9                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      10                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      11                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      12                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      13                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      14                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      15                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      16                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      17                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      18                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      19                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      20                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      21                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      22                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      23                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      24                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      25                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      26                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      27                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      28                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      29                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      30                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      31                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      32                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      33                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      34                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      35                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      36                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      37                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      38                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      39                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      40                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      41                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      42                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      43                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      44                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      45                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      46                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      47                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      48                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      49                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      50                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      51                 dispBeamColumn     {+0.0000 +1.0000 +0.0000}     {-1.0000 +0.0000 +0.0000}     {+0.0000 -0.0000 +1.0000}     { +Y -X +Z };   # Fiber-col
#      52                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      53                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      54                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      55                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      56                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      57                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      58                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      59                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      60                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      61                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      62                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      63                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      64                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      65                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      66                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      67                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      68                 dispBeamColumn     {+0.0000 -1.0000 +0.0000}     {+1.0000 +0.0000 -0.0000}     {+0.0000 +0.0000 +1.0000}     { -Y +X +Z };   # Fiber-col
#      69                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      70                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      71                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      72                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      73                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      74                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      75                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      76                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      77                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      78                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      79                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      80                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      81                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      82                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      83                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      84                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      85                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      86                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      87                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      88                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      89                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      90                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      91                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      92                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      93                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      94                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      95                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      96                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      97                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      98                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#      99                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     100                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     101                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     102                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     103                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     104                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     105                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     106                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     107                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     108                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     109                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     110                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     111                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     112                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
#     113                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-beam
