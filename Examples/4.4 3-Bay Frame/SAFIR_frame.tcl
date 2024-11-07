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

logFile "SAFIR_frame.log"
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

# Thermal_Displacement-Based_Beam-Column-Col (120)
# Thermal_Displacement-Based_Beam-Column-Beam-Mid (48)
# Thermal_Displacement-Based_Beam-Column-Beam-Joint (12)

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

node      1           18            0          9.6
node      2           18            0         9.28
node      3           18            0         8.96
node      4           18            0         8.64
node      5           18            0         8.32
node      6           18            0            8
node      7           18            0         7.68
node      8           18            0         7.36
node      9           18            0         7.04
node     10           18            0         6.72
node     11           18            0          6.4
node     12         17.4            0          6.4
node     13         16.8            0          6.4
node     14           18            0         6.08
node     15         16.2            0          6.4
node     16           18            0         5.76
node     17         15.6            0          6.4
node     18           18            0         5.44
node     19           15            0          6.4
node     20           18            0         5.12
node     21           18            0          4.8
node     22         14.4            0          6.4
node     23           18            0         4.48
node     24         13.8            0          6.4
node     25           18            0         4.16
node     26           18            0         3.84
node     27         13.2            0          6.4
node     28           12            0          9.6
node     29           12            0         9.28
node     30           12            0         8.96
node     31           12            0         8.64
node     32           18            0         3.52
node     33           12            0         8.32
node     34           12            0            8
node     35         12.6            0          6.4
node     36           12            0         7.68
node     37           18            0          3.2
node     38           12            0         7.36
node     39         17.4            0          3.2
node     40         16.8            0          3.2
node     41           12            0         7.04
node     42         16.2            0          3.2
node     43           12            0         6.72
node     44           18            0         2.88
node     45           12            0          6.4
node     46         15.6            0          3.2
node     47           12            0         6.08
node     48           18            0         2.56
node     49           15            0          3.2
node     50           12            0         5.76
node     51           12            0         5.44
node     52         11.4            0          6.4
node     53         14.4            0          3.2
node     54           18            0         2.24
node     55           12            0         5.12
node     56         13.8            0          3.2
node     57           18            0         1.92
node     58           12            0          4.8
node     59         10.8            0          6.4
node     60           12            0         4.48
node     61         13.2            0          3.2
node     62           18            0          1.6
node     63           12            0         4.16
node     64           12            0         3.84
node     65           18            0         1.28
node     66         12.6            0          3.2
node     67         10.2            0          6.4
node     68           12            0         3.52
node     69           18            0         0.96
node     70           12            0          3.2
node     71           18            0         0.64
node     72          9.6            0          6.4
node     73           12            0         2.88
node     74         11.4            0          3.2
node     75           12            0         2.56
node     76           18            0         0.32
node     77           12            0         2.24
node     78            9            0          6.4
node     79           18            0            0
node     80         10.8            0          3.2
node     81           12            0         1.92
node     82           12            0          1.6
node     83         10.2            0          3.2
node     84          8.4            0          6.4
node     85           12            0         1.28
node     86           12            0         0.96
node     87          9.6            0          3.2
node     88          7.8            0          6.4
node     89           12            0         0.64
node     90            9            0          3.2
node     91           12            0         0.32
node     92          7.2            0          6.4
node     93           12            0            0
node     94          8.4            0          3.2
node     95          6.6            0          6.4
node     96            6            0          9.6
node     97            6            0         9.28
node     98            6            0         8.96
node     99            6            0         8.64
node    100          7.8            0          3.2
node    101            6            0         8.32
node    102            6            0            8
node    103            6            0         7.68
node    104            6            0         7.36
node    105            6            0         7.04
node    106            6            0         6.72
node    107            6            0          6.4
node    108            6            0         6.08
node    109          7.2            0          3.2
node    110            6            0         5.76
node    111            6            0         5.44
node    112            6            0         5.12
node    113            6            0          4.8
node    114          5.4            0          6.4
node    115            6            0         4.48
node    116          6.6            0          3.2
node    117            6            0         4.16
node    118            6            0         3.84
node    119            6            0         3.52
node    120          4.8            0          6.4
node    121            6            0          3.2
node    122            6            0         2.88
node    123            6            0         2.56
node    124            6            0         2.24
node    125          5.4            0          3.2
node    126          4.2            0          6.4
node    127            6            0         1.92
node    128            6            0          1.6
node    129            6            0         1.28
node    130          4.8            0          3.2
node    131          3.6            0          6.4
node    132            6            0         0.96
node    133            6            0         0.64
node    134            6            0         0.32
node    135          4.2            0          3.2
node    136            3            0          6.4
node    137            6            0            0
node    138          3.6            0          3.2
node    139          2.4            0          6.4
node    140            3            0          3.2
node    141          1.8            0          6.4
node    142          2.4            0          3.2
node    143          1.2            0          6.4
node    144          1.8            0          3.2
node    145          0.6            0          6.4
node    146          1.2            0          3.2
node    147            0            0          9.6
node    148            0            0         9.28
node    149            0            0         8.96
node    150            0            0         8.64
node    151            0            0         8.32
node    152            0            0            8
node    153            0            0         7.68
node    154            0            0         7.36
node    155            0            0         7.04
node    156            0            0         6.72
node    157            0            0          6.4
node    158            0            0         6.08
node    159            0            0         5.76
node    160            0            0         5.44
node    161          0.6            0          3.2
node    162            0            0         5.12
node    163            0            0          4.8
node    164            0            0         4.48
node    165            0            0         4.16
node    166            0            0         3.84
node    167            0            0         3.52
node    168            0            0          3.2
node    169            0            0         2.88
node    170            0            0         2.56
node    171            0            0         2.24
node    172            0            0         1.92
node    173            0            0          1.6
node    174            0            0         1.28
node    175            0            0         0.96
node    176            0            0         0.64
node    177            0            0         0.32
node    178            0            0            0

# --------------------------------------------------------------------------------------------------------------
# R E S T R A I N T S
# --------------------------------------------------------------------------------------------------------------

# fix $NodeTag x-transl y-transl z-transl x-rot y-rot z-rot

fix      1   1   1   0   1   1   1
fix     28   1   1   0   1   1   1
fix     79   1   1   1   1   1   1
fix     93   1   1   1   1   1   1
fix     96   1   1   0   1   1   1
fix    137   1   1   1   1   1   1
fix    147   1   1   0   1   1   1
fix    178   1   1   1   1   1   1

fix    157   1   0 0 0 0 0
fix    168   1   0 0 0 0 0
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

### Section for element: Thermal_Displacement-Based_Beam-Column-Col
uniaxialMaterial DamagePlasticityConcreteECT 618 -30000 -0.0025 -200 -0.02 0.2 2900 2.72727e+06
uniaxialMaterial DPMsteelEC 611 500000 2e+08 0.05 0 1 0 1


section fiberSecThermal 667 -GJ 1e10  {


# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618     18     18  -0.160000  -0.160000   0.160000   0.160000

# Create the Cover fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618      3     24   0.160000  -0.200000   0.200000   0.200000
patch rect 618      3     24  -0.200000  -0.200000  -0.160000   0.200000
patch rect 618     18      3  -0.160000   0.160000   0.160000   0.200000
patch rect 618     18      3  -0.160000  -0.200000   0.160000  -0.160000

# Create the corner bars

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 611  2   0.00031420   0.160000   0.160000   0.160000  -0.160000
layer straight 611  2   0.00031420  -0.160000   0.160000  -0.160000  -0.160000

# Create the middle bars

# fiber $yLoc $zLoc $A $matTag
fiber   0.160000 0   0.00031420 611
fiber  -0.160000 0   0.00031420 611

# Create the middle Reinforcing Bars along local y axis

# fiber $yLoc $zLoc $A $matTag
fiber 0   0.160000     0.00031420 611
fiber 0  -0.160000     0.00031420 611
}
### Section for element: Thermal_Displacement-Based_Beam-Column-Beam-Mid


section fiberSecThermal 669 -GJ 1e10  {


# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618     13     15  -0.135000  -0.210000   0.135000   0.210000

# Create the Cover fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618      3     19   0.135000  -0.250000   0.175000   0.250000
patch rect 618      3     19  -0.175000  -0.250000  -0.135000   0.250000
patch rect 618     13      2  -0.135000   0.210000   0.135000   0.250000
patch rect 618     13      2  -0.135000  -0.250000   0.135000  -0.210000


# Create the Top bars (face on local z positive dir)

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 611   2    0.00031420  -0.135000   0.210000   0.135000   0.210000

# Create the Bottom bars (face on local z negative dir)

layer straight 611   3   0.00031420  -0.135000  -0.210000   0.135000  -0.210000
}
### Section for element: Thermal_Displacement-Based_Beam-Column-Beam-Joint


section fiberSecThermal 670 -GJ 1e10  {


# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618     13     15  -0.135000  -0.210000   0.135000   0.210000

# Create the Cover fibers

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618      3     19   0.135000  -0.250000   0.175000   0.250000
patch rect 618      3     19  -0.175000  -0.250000  -0.135000   0.250000
patch rect 618     13      2  -0.135000   0.210000   0.135000   0.250000
patch rect 618     13      2  -0.135000  -0.250000   0.135000  -0.210000


# Create the Top bars (face on local z positive dir)

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 611   3    0.00031420  -0.135000   0.210000   0.135000   0.210000

# Create the Bottom bars (face on local z negative dir)

layer straight 611   2   0.00031420  -0.135000  -0.210000   0.135000  -0.210000
}

# Displacement-Based Beam Column Element Definition

# element dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>

element dispBeamColumnThermal      1    178    177  3    667  5   -mass    1.256
element dispBeamColumnThermal      2    177    176  3    667  5   -mass    1.256
element dispBeamColumnThermal      3    176    175  3    667  5   -mass    1.256
element dispBeamColumnThermal      4    175    174  3    667  5   -mass    1.256
element dispBeamColumnThermal      5    174    173  3    667  5   -mass    1.256
element dispBeamColumnThermal      6    173    172  3    667  5   -mass    1.256
element dispBeamColumnThermal      7    172    171  3    667  5   -mass    1.256
element dispBeamColumnThermal      8    171    170  3    667  5   -mass    1.256
element dispBeamColumnThermal      9    170    169  3    667  5   -mass    1.256
element dispBeamColumnThermal     10    169    168  3    667  5   -mass    1.256
element dispBeamColumnThermal     11    168    167  3    667  5   -mass    1.256
element dispBeamColumnThermal     12    167    166  3    667  5   -mass    1.256
element dispBeamColumnThermal     13    166    165  3    667  5   -mass    1.256
element dispBeamColumnThermal     14    165    164  3    667  5   -mass    1.256
element dispBeamColumnThermal     15    164    163  3    667  5   -mass    1.256
element dispBeamColumnThermal     16    163    162  3    667  5   -mass    1.256
element dispBeamColumnThermal     17    162    160  3    667  5   -mass    1.256
element dispBeamColumnThermal     18    160    159  3    667  5   -mass    1.256
element dispBeamColumnThermal     19    159    158  3    667  5   -mass    1.256
element dispBeamColumnThermal     20    158    157  3    667  5   -mass    1.256
element dispBeamColumnThermal     21    157    156  3    667  5   -mass    1.256
element dispBeamColumnThermal     22    156    155  3    667  5   -mass    1.256
element dispBeamColumnThermal     23    155    154  3    667  5   -mass    1.256
element dispBeamColumnThermal     24    154    153  3    667  5   -mass    1.256
element dispBeamColumnThermal     25    153    152  3    667  5   -mass    1.256
element dispBeamColumnThermal     26    152    151  3    667  5   -mass    1.256
element dispBeamColumnThermal     27    151    150  3    667  5   -mass    1.256
element dispBeamColumnThermal     28    150    149  3    667  5   -mass    1.256
element dispBeamColumnThermal     29    149    148  3    667  5   -mass    1.256
element dispBeamColumnThermal     30    148    147  3    667  5   -mass    1.256
element dispBeamColumnThermal     31    168    161  3    670  6   -mass  1.37375
element dispBeamColumnThermal     32    161    146  3    669  6   -mass  1.37375
element dispBeamColumnThermal     33    146    144  3    669  6   -mass  1.37375
element dispBeamColumnThermal     34    144    142  3    669  6   -mass  1.37375
element dispBeamColumnThermal     35    142    140  3    669  6   -mass  1.37375
element dispBeamColumnThermal     36    140    138  3    669  6   -mass  1.37375
element dispBeamColumnThermal     37    138    135  3    669  6   -mass  1.37375
element dispBeamColumnThermal     38    135    130  3    669  6   -mass  1.37375
element dispBeamColumnThermal     39    130    125  3    669  6   -mass  1.37375
element dispBeamColumnThermal     40    125    121  3    670  6   -mass  1.37375
element dispBeamColumnThermal     41    121    116  3    670  6   -mass  1.37375
element dispBeamColumnThermal     42    116    109  3    669  6   -mass  1.37375
element dispBeamColumnThermal     43    109    100  3    669  6   -mass  1.37375
element dispBeamColumnThermal     44    100     94  3    669  6   -mass  1.37375
element dispBeamColumnThermal     45     94     90  3    669  6   -mass  1.37375
element dispBeamColumnThermal     46     90     87  3    669  6   -mass  1.37375
element dispBeamColumnThermal     47     87     83  3    669  6   -mass  1.37375
element dispBeamColumnThermal     48     83     80  3    669  6   -mass  1.37375
element dispBeamColumnThermal     49     80     74  3    669  6   -mass  1.37375
element dispBeamColumnThermal     50     74     70  3    670  6   -mass  1.37375
element dispBeamColumnThermal     51     70     66  3    670  6   -mass  1.37375
element dispBeamColumnThermal     52     66     61  3    669  6   -mass  1.37375
element dispBeamColumnThermal     53     61     56  3    669  6   -mass  1.37375
element dispBeamColumnThermal     54     56     53  3    669  6   -mass  1.37375
element dispBeamColumnThermal     55     53     49  3    669  6   -mass  1.37375
element dispBeamColumnThermal     56     49     46  3    669  6   -mass  1.37375
element dispBeamColumnThermal     57     46     42  3    669  6   -mass  1.37375
element dispBeamColumnThermal     58     42     40  3    669  6   -mass  1.37375
element dispBeamColumnThermal     59     40     39  3    669  6   -mass  1.37375
element dispBeamColumnThermal     60     39     37  3    670  6   -mass  1.37375
element dispBeamColumnThermal     61     37     32  3    667  5   -mass    1.256
element dispBeamColumnThermal     62     32     26  3    667  5   -mass    1.256
element dispBeamColumnThermal     63     26     25  3    667  5   -mass    1.256
element dispBeamColumnThermal     64     25     23  3    667  5   -mass    1.256
element dispBeamColumnThermal     65     23     21  3    667  5   -mass    1.256
element dispBeamColumnThermal     66     21     20  3    667  5   -mass    1.256
element dispBeamColumnThermal     67     20     18  3    667  5   -mass    1.256
element dispBeamColumnThermal     68     18     16  3    667  5   -mass    1.256
element dispBeamColumnThermal     69     16     14  3    667  5   -mass    1.256
element dispBeamColumnThermal     70     14     11  3    667  5   -mass    1.256
element dispBeamColumnThermal     71     11     10  3    667  5   -mass    1.256
element dispBeamColumnThermal     72     10      9  3    667  5   -mass    1.256
element dispBeamColumnThermal     73      9      8  3    667  5   -mass    1.256
element dispBeamColumnThermal     74      8      7  3    667  5   -mass    1.256
element dispBeamColumnThermal     75      7      6  3    667  5   -mass    1.256
element dispBeamColumnThermal     76      6      5  3    667  5   -mass    1.256
element dispBeamColumnThermal     77      5      4  3    667  5   -mass    1.256
element dispBeamColumnThermal     78      4      3  3    667  5   -mass    1.256
element dispBeamColumnThermal     79      3      2  3    667  5   -mass    1.256
element dispBeamColumnThermal     80      2      1  3    667  5   -mass    1.256
element dispBeamColumnThermal     81     28     29  3    667  5   -mass    1.256
element dispBeamColumnThermal     82     29     30  3    667  5   -mass    1.256
element dispBeamColumnThermal     83     30     31  3    667  5   -mass    1.256
element dispBeamColumnThermal     84     31     33  3    667  5   -mass    1.256
element dispBeamColumnThermal     85     33     34  3    667  5   -mass    1.256
element dispBeamColumnThermal     86     34     36  3    667  5   -mass    1.256
element dispBeamColumnThermal     87     36     38  3    667  5   -mass    1.256
element dispBeamColumnThermal     88     38     41  3    667  5   -mass    1.256
element dispBeamColumnThermal     89     41     43  3    667  5   -mass    1.256
element dispBeamColumnThermal     90     43     45  3    667  5   -mass    1.256
element dispBeamColumnThermal     91     45     47  3    667  5   -mass    1.256
element dispBeamColumnThermal     92     47     50  3    667  5   -mass    1.256
element dispBeamColumnThermal     93     50     51  3    667  5   -mass    1.256
element dispBeamColumnThermal     94     51     55  3    667  5   -mass    1.256
element dispBeamColumnThermal     95     55     58  3    667  5   -mass    1.256
element dispBeamColumnThermal     96     58     60  3    667  5   -mass    1.256
element dispBeamColumnThermal     97     60     63  3    667  5   -mass    1.256
element dispBeamColumnThermal     98     63     64  3    667  5   -mass    1.256
element dispBeamColumnThermal     99     64     68  3    667  5   -mass    1.256
element dispBeamColumnThermal    100     68     70  3    667  5   -mass    1.256
element dispBeamColumnThermal    101     70     73  3    667  5   -mass    1.256
element dispBeamColumnThermal    102     73     75  3    667  5   -mass    1.256
element dispBeamColumnThermal    103     75     77  3    667  5   -mass    1.256
element dispBeamColumnThermal    104     77     81  3    667  5   -mass    1.256
element dispBeamColumnThermal    105     81     82  3    667  5   -mass    1.256
element dispBeamColumnThermal    106     82     85  3    667  5   -mass    1.256
element dispBeamColumnThermal    107     85     86  3    667  5   -mass    1.256
element dispBeamColumnThermal    108     86     89  3    667  5   -mass    1.256
element dispBeamColumnThermal    109     89     91  3    667  5   -mass    1.256
element dispBeamColumnThermal    110     91     93  3    667  5   -mass    1.256
element dispBeamColumnThermal    111    137    134  3    667  5   -mass    1.256
element dispBeamColumnThermal    112    134    133  3    667  5   -mass    1.256
element dispBeamColumnThermal    113    133    132  3    667  5   -mass    1.256
element dispBeamColumnThermal    114    132    129  3    667  5   -mass    1.256
element dispBeamColumnThermal    115    129    128  3    667  5   -mass    1.256
element dispBeamColumnThermal    116    128    127  3    667  5   -mass    1.256
element dispBeamColumnThermal    117    127    124  3    667  5   -mass    1.256
element dispBeamColumnThermal    118    124    123  3    667  5   -mass    1.256
element dispBeamColumnThermal    119    123    122  3    667  5   -mass    1.256
element dispBeamColumnThermal    120    122    121  3    667  5   -mass    1.256
element dispBeamColumnThermal    121    121    119  3    667  5   -mass    1.256
element dispBeamColumnThermal    122    119    118  3    667  5   -mass    1.256
element dispBeamColumnThermal    123    118    117  3    667  5   -mass    1.256
element dispBeamColumnThermal    124    117    115  3    667  5   -mass    1.256
element dispBeamColumnThermal    125    115    113  3    667  5   -mass    1.256
element dispBeamColumnThermal    126    113    112  3    667  5   -mass    1.256
element dispBeamColumnThermal    127    112    111  3    667  5   -mass    1.256
element dispBeamColumnThermal    128    111    110  3    667  5   -mass    1.256
element dispBeamColumnThermal    129    110    108  3    667  5   -mass    1.256
element dispBeamColumnThermal    130    108    107  3    667  5   -mass    1.256
element dispBeamColumnThermal    131    107    106  3    667  5   -mass    1.256
element dispBeamColumnThermal    132    106    105  3    667  5   -mass    1.256
element dispBeamColumnThermal    133    105    104  3    667  5   -mass    1.256
element dispBeamColumnThermal    134    104    103  3    667  5   -mass    1.256
element dispBeamColumnThermal    135    103    102  3    667  5   -mass    1.256
element dispBeamColumnThermal    136    102    101  3    667  5   -mass    1.256
element dispBeamColumnThermal    137    101     99  3    667  5   -mass    1.256
element dispBeamColumnThermal    138     99     98  3    667  5   -mass    1.256
element dispBeamColumnThermal    139     98     97  3    667  5   -mass    1.256
element dispBeamColumnThermal    140     97     96  3    667  5   -mass    1.256
element dispBeamColumnThermal    141    157    145  3    670  6   -mass  1.37375
element dispBeamColumnThermal    142    145    143  3    669  6   -mass  1.37375
element dispBeamColumnThermal    143    143    141  3    669  6   -mass  1.37375
element dispBeamColumnThermal    144    141    139  3    669  6   -mass  1.37375
element dispBeamColumnThermal    145    139    136  3    669  6   -mass  1.37375
element dispBeamColumnThermal    146    136    131  3    669  6   -mass  1.37375
element dispBeamColumnThermal    147    131    126  3    669  6   -mass  1.37375
element dispBeamColumnThermal    148    126    120  3    669  6   -mass  1.37375
element dispBeamColumnThermal    149    120    114  3    669  6   -mass  1.37375
element dispBeamColumnThermal    150    114    107  3    670  6   -mass  1.37375
element dispBeamColumnThermal    151    107     95  3    670  6   -mass  1.37375
element dispBeamColumnThermal    152     95     92  3    669  6   -mass  1.37375
element dispBeamColumnThermal    153     92     88  3    669  6   -mass  1.37375
element dispBeamColumnThermal    154     88     84  3    669  6   -mass  1.37375
element dispBeamColumnThermal    155     84     78  3    669  6   -mass  1.37375
element dispBeamColumnThermal    156     78     72  3    669  6   -mass  1.37375
element dispBeamColumnThermal    157     72     67  3    669  6   -mass  1.37375
element dispBeamColumnThermal    158     67     59  3    669  6   -mass  1.37375
element dispBeamColumnThermal    159     59     52  3    669  6   -mass  1.37375
element dispBeamColumnThermal    160     52     45  3    670  6   -mass  1.37375
element dispBeamColumnThermal    161     45     35  3    670  6   -mass  1.37375
element dispBeamColumnThermal    162     35     27  3    669  6   -mass  1.37375
element dispBeamColumnThermal    163     27     24  3    669  6   -mass  1.37375
element dispBeamColumnThermal    164     24     22  3    669  6   -mass  1.37375
element dispBeamColumnThermal    165     22     19  3    669  6   -mass  1.37375
element dispBeamColumnThermal    166     19     17  3    669  6   -mass  1.37375
element dispBeamColumnThermal    167     17     15  3    669  6   -mass  1.37375
element dispBeamColumnThermal    168     15     13  3    669  6   -mass  1.37375
element dispBeamColumnThermal    169     13     12  3    669  6   -mass  1.37375
element dispBeamColumnThermal    170     12     11  3    670  6   -mass  1.37375
element dispBeamColumnThermal    171     37     44  3    667  5   -mass    1.256
element dispBeamColumnThermal    172     44     48  3    667  5   -mass    1.256
element dispBeamColumnThermal    173     48     54  3    667  5   -mass    1.256
element dispBeamColumnThermal    174     54     57  3    667  5   -mass    1.256
element dispBeamColumnThermal    175     57     62  3    667  5   -mass    1.256
element dispBeamColumnThermal    176     62     65  3    667  5   -mass    1.256
element dispBeamColumnThermal    177     65     69  3    667  5   -mass    1.256
element dispBeamColumnThermal    178     69     71  3    667  5   -mass    1.256
element dispBeamColumnThermal    179     71     76  3    667  5   -mass    1.256
element dispBeamColumnThermal    180     76     79  3    667  5   -mass    1.256

# --------------------------------------------------------------------------------------------------------------
#
# D O M A I N  C O M M O N S
#
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# R E C O R D E R S
# --------------------------------------------------------------------------------------------------------------

recorder Node -file Node_displacements.out -time -nodeRange 1 178 -dof 1 2 3 disp
recorder Node -file Node_rotations.out -time -nodeRange 1 178 -dof 4 5 6 disp
recorder Node -file Node_forceReactions.out -time -nodeRange 1 178 -dof 1 2 3 reaction
recorder Node -file Node_momentReactions.out -time -nodeRange 1 178 -dof 4 5 6 reaction
recorder Element -file DispBeamColumn_localForce.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 localForce
recorder Element -file DispBeamColumn_basicDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 basicDeformation
recorder Element -file DispBeamColumn_plasticDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 plasticDeformation
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
    load      1        0        0     -880        0        0        0
    load     28        0        0    -1760        0        0        0
    load     96        0        0    -1760        0        0        0
    load    147        0        0     -880        0        0        0
    eleLoad -ele     31 -type -beamUniform        0      -42        0
    eleLoad -ele     32 -type -beamUniform        0      -42        0
    eleLoad -ele     33 -type -beamUniform        0      -42        0
    eleLoad -ele     34 -type -beamUniform        0      -42        0
    eleLoad -ele     35 -type -beamUniform        0      -42        0
    eleLoad -ele     36 -type -beamUniform        0      -42        0
    eleLoad -ele     37 -type -beamUniform        0      -42        0
    eleLoad -ele     38 -type -beamUniform        0      -42        0
    eleLoad -ele     39 -type -beamUniform        0      -42        0
    eleLoad -ele     40 -type -beamUniform        0      -42        0
    eleLoad -ele     41 -type -beamUniform        0      -42        0
    eleLoad -ele     42 -type -beamUniform        0      -42        0
    eleLoad -ele     43 -type -beamUniform        0      -42        0
    eleLoad -ele     44 -type -beamUniform        0      -42        0
    eleLoad -ele     45 -type -beamUniform        0      -42        0
    eleLoad -ele     46 -type -beamUniform        0      -42        0
    eleLoad -ele     47 -type -beamUniform        0      -42        0
    eleLoad -ele     48 -type -beamUniform        0      -42        0
    eleLoad -ele     49 -type -beamUniform        0      -42        0
    eleLoad -ele     50 -type -beamUniform        0      -42        0
    eleLoad -ele     51 -type -beamUniform        0      -42        0
    eleLoad -ele     52 -type -beamUniform        0      -42        0
    eleLoad -ele     53 -type -beamUniform        0      -42        0
    eleLoad -ele     54 -type -beamUniform        0      -42        0
    eleLoad -ele     55 -type -beamUniform        0      -42        0
    eleLoad -ele     56 -type -beamUniform        0      -42        0
    eleLoad -ele     57 -type -beamUniform        0      -42        0
    eleLoad -ele     58 -type -beamUniform        0      -42        0
    eleLoad -ele     59 -type -beamUniform        0      -42        0
    eleLoad -ele     60 -type -beamUniform        0      -42        0
    eleLoad -ele    141 -type -beamUniform        0      -42        0
    eleLoad -ele    142 -type -beamUniform        0      -42        0
    eleLoad -ele    143 -type -beamUniform        0      -42        0
    eleLoad -ele    144 -type -beamUniform        0      -42        0
    eleLoad -ele    145 -type -beamUniform        0      -42        0
    eleLoad -ele    146 -type -beamUniform        0      -42        0
    eleLoad -ele    147 -type -beamUniform        0      -42        0
    eleLoad -ele    148 -type -beamUniform        0      -42        0
    eleLoad -ele    149 -type -beamUniform        0      -42        0
    eleLoad -ele    150 -type -beamUniform        0      -42        0
    eleLoad -ele    151 -type -beamUniform        0      -42        0
    eleLoad -ele    152 -type -beamUniform        0      -42        0
    eleLoad -ele    153 -type -beamUniform        0      -42        0
    eleLoad -ele    154 -type -beamUniform        0      -42        0
    eleLoad -ele    155 -type -beamUniform        0      -42        0
    eleLoad -ele    156 -type -beamUniform        0      -42        0
    eleLoad -ele    157 -type -beamUniform        0      -42        0
    eleLoad -ele    158 -type -beamUniform        0      -42        0
    eleLoad -ele    159 -type -beamUniform        0      -42        0
    eleLoad -ele    160 -type -beamUniform        0      -42        0
    eleLoad -ele    161 -type -beamUniform        0      -42        0
    eleLoad -ele    162 -type -beamUniform        0      -42        0
    eleLoad -ele    163 -type -beamUniform        0      -42        0
    eleLoad -ele    164 -type -beamUniform        0      -42        0
    eleLoad -ele    165 -type -beamUniform        0      -42        0
    eleLoad -ele    166 -type -beamUniform        0      -42        0
    eleLoad -ele    167 -type -beamUniform        0      -42        0
    eleLoad -ele    168 -type -beamUniform        0      -42        0
    eleLoad -ele    169 -type -beamUniform        0      -42        0
    eleLoad -ele    170 -type -beamUniform        0      -42        0

}

# recording the initial status

record

# Analysis options

system BandGeneral
numberer RCM
constraints Transformation
integrator LoadControl 0.01
test RelativeEnergyIncr 1e-05 50 2
algorithm Newton
analysis Static
set Lincr 0.01
set Nsteps 100
set committedSteps 1
set LoadCounter 0


set strIni {}
variable testTypeStatic RelativeEnergyIncr
variable TolStatic 1e-05
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
# Loads - Plain Pattern
pattern Plain 200 Linear {
    eleLoad -ele      1 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      2 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      3 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      4 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      5 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      6 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      7 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      8 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele      9 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     10 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     11 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     12 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     13 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     14 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     15 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     16 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     17 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     18 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     19 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     20 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     31 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     32 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     33 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     34 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     35 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     36 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     37 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     38 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     39 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     40 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     41 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     42 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     43 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     44 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     45 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     46 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     47 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     48 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     49 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     50 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     51 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     52 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     53 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     54 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     55 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     56 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     57 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     58 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     59 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     60 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25     
    eleLoad -ele     61 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     62 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     63 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     64 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     65 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     66 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     67 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     68 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     69 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     70 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     91 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     92 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     93 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     94 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     95 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     96 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     97 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     98 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele     99 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    100 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    101 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    102 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    103 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    104 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    105 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    106 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    107 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    108 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    109 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    110 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    111 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    112 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    113 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    114 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    115 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    116 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    117 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    118 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    119 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    120 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    121 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    122 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    123 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    124 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    125 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    126 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    127 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    128 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    129 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    130 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    141 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    142 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    143 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    144 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    145 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    146 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    147 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    148 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    149 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    150 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    151 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    152 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    153 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    154 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    155 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    156 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    157 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    158 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    159 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    160 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    161 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    162 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    163 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    164 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    165 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    166 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    167 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    168 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    169 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    170 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation400   -0.175    0.175    -0.25     0.25      
    eleLoad -ele    171 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    172 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    173 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    174 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    175 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    176 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    177 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    178 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    179 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24
    eleLoad -ele    180 -type -beamThermal -z   -source "ColumnRC1.dat"  -genInterpolation     -0.2      0.2     -0.2      0.2       24       24

}

# recording the initial status

record

# Analysis options

system UmfPack
numberer RCM
constraints Transformation
integrator LoadControl 30
test RelativeEnergyIncr 0.0001 50 2
algorithm Newton
analysis Static
set Lincr 30
set Nsteps 420
set committedSteps 1
set LoadCounter 0


set strIni {}
variable testTypeStatic RelativeEnergyIncr
variable TolStatic 0.0001
variable maxNumIterStatic 50
variable algorithmTypeStatic Newton

for {set i 1} { $i <= $Nsteps } {incr i 1} {
    set t [format "%7.5f" [expr [getTime] + $Lincr]]
    puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
                puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
                puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
                puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
                puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
                puts -nonewline "(2) $algorithmTypeStatic$strIni LF $t "
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
# 178

# Elements 1D
# 180

# Elements 2D
# 0

# Elements 3D
# 0

# DispBeamColumn
# 180
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# F R A M E   L O C A L   A X E S   O R I E N T A T I O N
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#      ID                           Type                       Local-x                       Local-y                       Local-z          Literal      Material / Section
#
#       1                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       2                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       3                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       4                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       5                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       6                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       7                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       8                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#       9                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      10                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      11                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      12                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      13                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      14                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      15                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      16                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      17                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      18                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      19                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      20                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      21                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      22                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      23                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      24                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      25                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      26                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      27                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      28                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      29                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      30                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      31                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      32                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      33                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      34                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      35                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      36                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      37                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      38                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      39                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      40                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      41                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      42                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      43                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      44                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      45                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      46                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      47                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      48                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      49                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      50                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      51                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      52                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      53                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      54                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      55                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      56                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      57                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      58                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      59                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#      60                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#      61                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      62                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      63                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      64                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      65                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      66                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      67                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      68                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      69                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      70                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      71                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      72                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      73                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      74                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      75                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      76                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      77                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      78                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      79                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      80                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#      81                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      82                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      83                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      84                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      85                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      86                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      87                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      88                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      89                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      90                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      91                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      92                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      93                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      94                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      95                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      96                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      97                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      98                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#      99                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     100                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     101                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     102                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     103                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     104                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     105                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     106                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     107                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     108                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     109                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     110                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     111                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     112                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     113                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     114                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     115                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     116                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     117                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     118                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     119                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     120                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     121                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     122                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     123                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     124                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     125                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     126                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     127                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     128                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     129                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     130                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     131                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     132                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     133                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     134                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     135                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     136                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     137                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     138                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     139                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     140                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber-Col
#     141                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     142                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     143                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     144                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     145                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     146                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     147                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     148                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     149                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     150                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     151                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     152                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     153                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     154                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     155                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     156                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     157                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     158                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     159                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     160                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     161                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     162                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     163                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     164                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     165                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     166                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     167                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     168                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     169                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Mid
#     170                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber-Beam-Joint
#     171                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     172                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     173                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     174                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     175                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     176                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     177                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     178                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     179                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
#     180                 dispBeamColumn     {+0.0000 +0.0000 -1.0000}     {-0.0000 -1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { -Z -Y -X };   # Fiber-Col
