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

logFile "C1Gernay2022.log"
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

# Thermal_Displacement-Based_Beam-Column (90)

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

node      1            0            0            3
node      2            0            0      2.96667
node      3            0            0      2.93333
node      4            0            0          2.9
node      5            0            0      2.86667
node      6            0            0      2.83333
node      7            0            0          2.8
node      8            0            0      2.76667
node      9            0            0      2.73333
node     10            0            0          2.7
node     11            0            0      2.66667
node     12            0            0      2.63333
node     13            0            0          2.6
node     14            0            0      2.56667
node     15            0            0      2.53333
node     16            0            0          2.5
node     17            0            0      2.46667
node     18            0            0      2.43333
node     19            0            0          2.4
node     20            0            0      2.36667
node     21            0            0      2.33333
node     22            0            0          2.3
node     23            0            0      2.26667
node     24            0            0      2.23333
node     25            0            0          2.2
node     26            0            0      2.16667
node     27            0            0      2.13333
node     28            0            0          2.1
node     29            0            0      2.06667
node     30            0            0      2.03333
node     31            0            0            2
node     32            0            0      1.96667
node     33            0            0      1.93333
node     34            0            0          1.9
node     35            0            0      1.86667
node     36            0            0      1.83333
node     37            0            0          1.8
node     38            0            0      1.76667
node     39            0            0      1.73333
node     40            0            0          1.7
node     41            0            0      1.66667
node     42            0            0      1.63333
node     43            0            0          1.6
node     44            0            0      1.56667
node     45            0            0      1.53333
node     46            0            0          1.5
node     47            0            0      1.46667
node     48            0            0      1.43333
node     49            0            0          1.4
node     50            0            0      1.36667
node     51            0            0      1.33333
node     52            0            0          1.3
node     53            0            0      1.26667
node     54            0            0      1.23333
node     55            0            0          1.2
node     56            0            0      1.16667
node     57            0            0      1.13333
node     58            0            0          1.1
node     59            0            0      1.06667
node     60            0            0      1.03333
node     61            0            0            1
node     62            0            0     0.966667
node     63            0            0     0.933333
node     64            0            0          0.9
node     65            0            0     0.866667
node     66            0            0     0.833333
node     67            0            0          0.8
node     68            0            0     0.766667
node     69            0            0     0.733333
node     70            0            0          0.7
node     71            0            0     0.666667
node     72            0            0     0.633333
node     73            0            0          0.6
node     74            0            0     0.566667
node     75            0            0     0.533333
node     76            0            0          0.5
node     77            0            0     0.466667
node     78            0            0     0.433333
node     79            0            0          0.4
node     80            0            0     0.366667
node     81            0            0     0.333333
node     82            0            0          0.3
node     83            0            0     0.266667
node     84            0            0     0.233333
node     85            0            0          0.2
node     86            0            0     0.166667
node     87            0            0     0.133333
node     88            0            0          0.1
node     89            0            0    0.0666667
node     90            0            0    0.0333333
node     91            0            0            0

# --------------------------------------------------------------------------------------------------------------
# R E S T R A I N T S
# --------------------------------------------------------------------------------------------------------------

# fix $NodeTag x-transl y-transl z-transl x-rot y-rot z-rot
fix      1   1   1   0   1   0   1
fix     91   1   1   1   0   0   0

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

### Section for element: Thermal_Displacement-Based_Beam-Column
uniaxialMaterial DamagePlasticityConcreteECT 667 -34000 -0.0025 -3400 -0.02 0.1 3400 3.42727e+06
uniaxialMaterial DPMsteelEC 668 500000 2e+08 0.01 0 1 0 1


section fiberSecThermal 642 -GJ 1e10  {

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 667     30     30  -0.15000  -0.15000   0.15000   0.15000

# Create the corner bars

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 668  2   0.00015390   0.117000   0.117000   0.117000  -0.117000
layer straight 668  2   0.00015390  -0.117000   0.117000  -0.117000  -0.117000

# Create the middle bars

# fiber $yLoc $zLoc $A $matTag
fiber   0.119000 0   0.00007854 668
fiber  -0.119000 0   0.00007854 668

# Create the middle Reinforcing Bars along local y axis

# fiber $yLoc $zLoc $A $matTag
fiber 0   0.119000     0.00007854 668
fiber 0  -0.119000     0.00007854 668
}

# Displacement-Based Beam Column Element Definition

# element dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>

element dispBeamColumnThermal      1     91     90  3    642  5   -mass    7.065
element dispBeamColumnThermal      2     90     89  3    642  5   -mass    7.065
element dispBeamColumnThermal      3     89     88  3    642  5   -mass    7.065
element dispBeamColumnThermal      4     88     87  3    642  5   -mass    7.065
element dispBeamColumnThermal      5     87     86  3    642  5   -mass    7.065
element dispBeamColumnThermal      6     86     85  3    642  5   -mass    7.065
element dispBeamColumnThermal      7     85     84  3    642  5   -mass    7.065
element dispBeamColumnThermal      8     84     83  3    642  5   -mass    7.065
element dispBeamColumnThermal      9     83     82  3    642  5   -mass    7.065
element dispBeamColumnThermal     10     82     81  3    642  5   -mass    7.065
element dispBeamColumnThermal     11     81     80  3    642  5   -mass    7.065
element dispBeamColumnThermal     12     80     79  3    642  5   -mass    7.065
element dispBeamColumnThermal     13     79     78  3    642  5   -mass    7.065
element dispBeamColumnThermal     14     78     77  3    642  5   -mass    7.065
element dispBeamColumnThermal     15     77     76  3    642  5   -mass    7.065
element dispBeamColumnThermal     16     76     75  3    642  5   -mass    7.065
element dispBeamColumnThermal     17     75     74  3    642  5   -mass    7.065
element dispBeamColumnThermal     18     74     73  3    642  5   -mass    7.065
element dispBeamColumnThermal     19     73     72  3    642  5   -mass    7.065
element dispBeamColumnThermal     20     72     71  3    642  5   -mass    7.065
element dispBeamColumnThermal     21     71     70  3    642  5   -mass    7.065
element dispBeamColumnThermal     22     70     69  3    642  5   -mass    7.065
element dispBeamColumnThermal     23     69     68  3    642  5   -mass    7.065
element dispBeamColumnThermal     24     68     67  3    642  5   -mass    7.065
element dispBeamColumnThermal     25     67     66  3    642  5   -mass    7.065
element dispBeamColumnThermal     26     66     65  3    642  5   -mass    7.065
element dispBeamColumnThermal     27     65     64  3    642  5   -mass    7.065
element dispBeamColumnThermal     28     64     63  3    642  5   -mass    7.065
element dispBeamColumnThermal     29     63     62  3    642  5   -mass    7.065
element dispBeamColumnThermal     30     62     61  3    642  5   -mass    7.065
element dispBeamColumnThermal     31     61     60  3    642  5   -mass    7.065
element dispBeamColumnThermal     32     60     59  3    642  5   -mass    7.065
element dispBeamColumnThermal     33     59     58  3    642  5   -mass    7.065
element dispBeamColumnThermal     34     58     57  3    642  5   -mass    7.065
element dispBeamColumnThermal     35     57     56  3    642  5   -mass    7.065
element dispBeamColumnThermal     36     56     55  3    642  5   -mass    7.065
element dispBeamColumnThermal     37     55     54  3    642  5   -mass    7.065
element dispBeamColumnThermal     38     54     53  3    642  5   -mass    7.065
element dispBeamColumnThermal     39     53     52  3    642  5   -mass    7.065
element dispBeamColumnThermal     40     52     51  3    642  5   -mass    7.065
element dispBeamColumnThermal     41     51     50  3    642  5   -mass    7.065
element dispBeamColumnThermal     42     50     49  3    642  5   -mass    7.065
element dispBeamColumnThermal     43     49     48  3    642  5   -mass    7.065
element dispBeamColumnThermal     44     48     47  3    642  5   -mass    7.065
element dispBeamColumnThermal     45     47     46  3    642  5   -mass    7.065
element dispBeamColumnThermal     46     46     45  3    642  5   -mass    7.065
element dispBeamColumnThermal     47     45     44  3    642  5   -mass    7.065
element dispBeamColumnThermal     48     44     43  3    642  5   -mass    7.065
element dispBeamColumnThermal     49     43     42  3    642  5   -mass    7.065
element dispBeamColumnThermal     50     42     41  3    642  5   -mass    7.065
element dispBeamColumnThermal     51     41     40  3    642  5   -mass    7.065
element dispBeamColumnThermal     52     40     39  3    642  5   -mass    7.065
element dispBeamColumnThermal     53     39     38  3    642  5   -mass    7.065
element dispBeamColumnThermal     54     38     37  3    642  5   -mass    7.065
element dispBeamColumnThermal     55     37     36  3    642  5   -mass    7.065
element dispBeamColumnThermal     56     36     35  3    642  5   -mass    7.065
element dispBeamColumnThermal     57     35     34  3    642  5   -mass    7.065
element dispBeamColumnThermal     58     34     33  3    642  5   -mass    7.065
element dispBeamColumnThermal     59     33     32  3    642  5   -mass    7.065
element dispBeamColumnThermal     60     32     31  3    642  5   -mass    7.065
element dispBeamColumnThermal     61     31     30  3    642  5   -mass    7.065
element dispBeamColumnThermal     62     30     29  3    642  5   -mass    7.065
element dispBeamColumnThermal     63     29     28  3    642  5   -mass    7.065
element dispBeamColumnThermal     64     28     27  3    642  5   -mass    7.065
element dispBeamColumnThermal     65     27     26  3    642  5   -mass    7.065
element dispBeamColumnThermal     66     26     25  3    642  5   -mass    7.065
element dispBeamColumnThermal     67     25     24  3    642  5   -mass    7.065
element dispBeamColumnThermal     68     24     23  3    642  5   -mass    7.065
element dispBeamColumnThermal     69     23     22  3    642  5   -mass    7.065
element dispBeamColumnThermal     70     22     21  3    642  5   -mass    7.065
element dispBeamColumnThermal     71     21     20  3    642  5   -mass    7.065
element dispBeamColumnThermal     72     20     19  3    642  5   -mass    7.065
element dispBeamColumnThermal     73     19     18  3    642  5   -mass    7.065
element dispBeamColumnThermal     74     18     17  3    642  5   -mass    7.065
element dispBeamColumnThermal     75     17     16  3    642  5   -mass    7.065
element dispBeamColumnThermal     76     16     15  3    642  5   -mass    7.065
element dispBeamColumnThermal     77     15     14  3    642  5   -mass    7.065
element dispBeamColumnThermal     78     14     13  3    642  5   -mass    7.065
element dispBeamColumnThermal     79     13     12  3    642  5   -mass    7.065
element dispBeamColumnThermal     80     12     11  3    642  5   -mass    7.065
element dispBeamColumnThermal     81     11     10  3    642  5   -mass    7.065
element dispBeamColumnThermal     82     10      9  3    642  5   -mass    7.065
element dispBeamColumnThermal     83      9      8  3    642  5   -mass    7.065
element dispBeamColumnThermal     84      8      7  3    642  5   -mass    7.065
element dispBeamColumnThermal     85      7      6  3    642  5   -mass    7.065
element dispBeamColumnThermal     86      6      5  3    642  5   -mass    7.065
element dispBeamColumnThermal     87      5      4  3    642  5   -mass    7.065
element dispBeamColumnThermal     88      4      3  3    642  5   -mass    7.065
element dispBeamColumnThermal     89      3      2  3    642  5   -mass    7.065
element dispBeamColumnThermal     90      2      1  3    642  5   -mass    7.065

# --------------------------------------------------------------------------------------------------------------
#
# D O M A I N  C O M M O N S
#
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# R E C O R D E R S
# --------------------------------------------------------------------------------------------------------------

recorder Node -file Node_displacements.out -time -nodeRange 1 91 -dof 1 2 3 disp
recorder Node -file Node_rotations.out -time -nodeRange 1 91 -dof 4 5 6 disp
recorder Node -file Node_forceReactions.out -time -nodeRange 1 91 -dof 1 2 3 reaction
recorder Node -file Node_momentReactions.out -time -nodeRange 1 91 -dof 4 5 6 reaction
recorder Element -file DispBeamColumn_localForce.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 localForce
recorder Element -file DispBeamColumn_basicDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 basicDeformation
recorder Element -file DispBeamColumn_plasticDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 plasticDeformation
puts ""
puts " __   __       __          __                   _       "
puts "/ _ .|  \\ _|_ /  \\ _  _ _ (_  _ _ _  | _ |_ _ _(_ _  _ _"
puts "\\__)||__/  |  \\__/|_)(-| )__)(-(-_)  || )|_(-| | (_|(_(-"
puts "                  |                                     "
puts "                                                  v2.8.0"
puts "Analysis summary"
puts ""
puts "Interval 1 : Static - [expr int(1+100)] steps"
puts "Interval 2 : Transient - [expr int(1.0 + 10/0.005)] steps x 0.005 s"
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
    load      1        0        0    -1009      0     20     0
	#load      46       -29        0     0         0     0        0
	
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
# Loads - Plain Pattern
pattern Plain 200 Linear {
    eleLoad -ele      1 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      2 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      3 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      4 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      5 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      6 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      7 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      8 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele      9 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     10 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     11 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     12 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     13 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     14 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     15 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     16 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     17 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     18 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     19 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     20 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15      0.15       30  30
    eleLoad -ele     21 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     22 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     23 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     24 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     25 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     26 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     27 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     28 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     29 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     30 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     31 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     32 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     33 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     34 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     35 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     36 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     37 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     38 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     39 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     40 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     41 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     42 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     43 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     44 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     45 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     46 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     47 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     48 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     49 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     50 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     51 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     52 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     53 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     54 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     55 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     56 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     57 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     58 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     59 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     60 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     61 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     62 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     63 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     64 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     65 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     66 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     67 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     68 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     69 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     70 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     71 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     72 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     73 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     74 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     75 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     76 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     77 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     78 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     79 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     80 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     81 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     82 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     83 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     84 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     85 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     86 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     87 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     88 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     89 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30
    eleLoad -ele     90 -type -beamThermal -z   -source "BeamColumnRCT4.dat"  -genInterpolation     -0.15      0.15     -0.15     0.15        30  30

}
# recording the initial status

record

# Analysis options

system UmfPack
numberer RCM
constraints Transformation
integrator Newmark 0.45 0.8
test RelativeEnergyIncr 0.0001 500 2
algorithm KrylovNewton
analysis Transient


# Loads - Uniform Excitation
set iGMdirection 1
set DtAnalysis 300
set TmaxAnalysis 9000

set committedSteps 1
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)]

set strIni {}
variable testTypeDynamic RelativeEnergyIncr
variable TolDynamic 0.0001
variable maxNumIterDynamic 500
variable algorithmTypeDynamic NewtonLineSearch

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
#remove loadPattern 200

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

puts "Running interval 3\n"
# Loads - Plain Pattern
pattern Plain 700 Linear {
    load      1        0        0    -8        0     0        0
	
}
# recording the initial status

record

# Analysis options

system UmfPack
numberer RCM
constraints Transformation
integrator Newmark 0.5 0.25
test RelativeEnergyIncr 0.0001 500 2
algorithm KrylovNewton
analysis Transient


# Loads - Uniform Excitation
set iGMdirection 1
set DtAnalysis 0.01
set TmaxAnalysis 1

set committedSteps 1
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)]

set strIni {}
variable testTypeDynamic RelativeEnergyIncr
variable TolDynamic 0.0001
variable maxNumIterDynamic 500
variable algorithmTypeDynamic NewtonLineSearch

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
#remove loadPattern 200

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






# Number of nodes
# 91

# Elements 1D
# 90

# Elements 2D
# 0

# Elements 3D
# 0

# DispBeamColumn
# 90
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# F R A M E   L O C A L   A X E S   O R I E N T A T I O N
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#      ID                           Type                       Local-x                       Local-y                       Local-z          Literal      Material / Section
#
#       1                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       2                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       3                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       4                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       5                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       6                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       7                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       8                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#       9                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      10                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      11                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      12                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      13                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      14                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      15                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      16                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      17                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      18                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      19                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      20                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      21                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      22                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      23                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      24                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      25                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      26                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      27                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      28                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      29                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      30                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      31                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      32                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      33                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      34                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      35                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      36                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      37                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      38                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      39                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      40                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      41                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      42                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      43                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      44                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      45                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      46                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      47                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      48                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      49                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      50                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      51                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      52                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      53                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      54                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      55                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      56                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      57                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      58                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      59                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      60                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      61                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      62                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      63                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      64                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      65                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      66                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      67                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      68                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      69                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      70                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      71                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      72                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      73                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      74                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      75                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      76                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      77                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      78                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      79                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      80                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      81                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      82                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      83                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      84                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      85                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      86                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      87                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      88                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      89                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
#      90                 dispBeamColumn     {+0.0000 +0.0000 +1.0000}     {+0.0000 +1.0000 -0.0000}     {-1.0000 +0.0000 +0.0000}     { +Z +Y -X };   # Fiber
