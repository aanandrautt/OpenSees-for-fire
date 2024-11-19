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

logFile "Beam.log"
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

# Thermal_Displacement-Based_Beam-Column (14)

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

node      1          6.1            0            0
node      2      5.59167            0            0
node      3      5.08333            0            0
node      4        4.575            0            0
node      5      4.06667            0            0
node      6      3.55833            0            0
node      7         3.05            0            0
node      9      2.54167            0            0
node     10      2.03333            0            0
node     11        1.525            0            0
node     13      1.01667            0            0
node     14     0.508333            0            0
node     15            0            0            0


mass 4 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 5 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 6 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 7 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 9 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 10 1e-9 1e-9 -25 1e-9 1e-9 1e-9
mass 11 1e-9 1e-9 -25 1e-9 1e-9 1e-9


# --------------------------------------------------------------------------------------------------------------
# R E S T R A I N T S
# --------------------------------------------------------------------------------------------------------------

# fix $NodeTag x-transl y-transl z-transl x-rot y-rot z-rot

fix      1   0   1   1   0   0   0
fix     15   1   1   1   0   0   0

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
uniaxialMaterial DamagePlasticityConcreteECT 618 -15000 -0.0025 -11000 -0.02 0.2 2200 2e+06
uniaxialMaterial DPMsteelEC 611 310000 2e+08 0.1 0 1 0 1


section fiberSecThermal 642 -GJ 1e10  {

# patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
patch rect 618      16     24  -0.10000  -0.30000   0.10000   0.30000

# Create the Top bars (face on local z positive dir)

# layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
layer straight 611   2    0.00011310  -0.060000   0.260000   0.060000   0.260000

# Create the Bottom bars (face on local z negative dir)

layer straight 611   3   0.00038010  -0.060000  -0.260000   0.060000  -0.260000
}

# Displacement-Based Beam Column Element Definition

# element dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>

element dispBeamColumnThermal      1     15     14  3    642  6   -mass    0.942
element dispBeamColumnThermal      2     14     13  3    642  6   -mass    0.942
element dispBeamColumnThermal      3     13     11  3    642  6   -mass    0.942
element dispBeamColumnThermal      4     11     10  3    642  6   -mass    0.942
element dispBeamColumnThermal      5     10      9  3    642  6   -mass    0.942
element dispBeamColumnThermal      6      9      7  3    642  6   -mass    0.942
element dispBeamColumnThermal      7      7      6  3    642  6   -mass    0.942
element dispBeamColumnThermal      8      6      5  3    642  6   -mass    0.942
element dispBeamColumnThermal      9      5      4  3    642  6   -mass    0.942
element dispBeamColumnThermal     10      4      3  3    642  6   -mass    0.942
element dispBeamColumnThermal     11      3      2  3    642  6   -mass    0.942
element dispBeamColumnThermal     12      2      1  3    642  6   -mass    0.942

# --------------------------------------------------------------------------------------------------------------
#
# D O M A I N  C O M M O N S
#
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# R E C O R D E R S
# --------------------------------------------------------------------------------------------------------------

recorder Node -file Node_displacements.out -time -nodeRange 1 15 -dof 1 2 3 disp
recorder Node -file Node_rotations.out -time -nodeRange 1 15 -dof 4 5 6 disp
recorder Node -file Node_forceReactions.out -time -nodeRange 1 15 -dof 1 2 3 reaction
recorder Node -file Node_momentReactions.out -time -nodeRange 1 15 -dof 4 5 6 reaction
recorder Element -file DispBeamColumn_localForce.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 localForce
recorder Element -file DispBeamColumn_basicDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 basicDeformation
recorder Element -file DispBeamColumn_plasticDeformation.out -time -ele 1 2 3 4 5 6 7 8 9 10 11 12 13 14 plasticDeformation

puts ""
puts " __   __       __          __                   _       "
puts "/ _ .|  \\ _|_ /  \\ _  _ _ (_  _ _ _  | _ |_ _ _(_ _  _ _"
puts "\\__)||__/  |  \\__/|_)(-| )__)(-(-_)  || )|_(-| | (_|(_(-"
puts "                  |                                     "
puts "                                                  v2.8.0"
puts "Analysis summary"
puts ""
puts "Interval 1 : Static - [expr int(1+50)] steps"
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
    load      4        0       0   -32.5     0        0        0
    load     11        0       0   -32.5     0        0        0

}

# recording the initial status

record

# Analysis options

system BandGeneral
numberer RCM
constraints Transformation
integrator LoadControl 0.02
test RelativeEnergyIncr 0.001 50 2
algorithm Newton
analysis Static
set Lincr 0.02
set Nsteps 50
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
    eleLoad -ele      1 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      2 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      3 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      4 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      5 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      6 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      7 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      8 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele      9 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele     10 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele     11 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
    eleLoad -ele     12 -type -beamThermal -z   -source "BeamRC1.dat"  -genInterpolation     -0.1      0.1     -0.3      0.3       16  24
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
set DtAnalysis 30
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
wipe;
# --------------------------------------------------------------------------------------------------------------
#
# M E T A D A T A
#
# --------------------------------------------------------------------------------------------------------------

# Number of nodes
# 15

# Elements 1D
# 14

# Elements 2D
# 0

# Elements 3D
# 0

# DispBeamColumn
# 14
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# F R A M E   L O C A L   A X E S   O R I E N T A T I O N
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#      ID                           Type                       Local-x                       Local-y                       Local-z          Literal      Material / Section
#
#       1                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       2                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       3                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       4                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       5                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       6                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       7                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       8                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#       9                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#      10                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#      11                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#      12                 dispBeamColumn     {+1.0000 +0.0000 +0.0000}     {+0.0000 +1.0000 +0.0000}     {+0.0000 +0.0000 +1.0000}     { +X +Y +Z };   # Fiber
#      13                 dispBeamColumn     {+0.9487 +0.0000 +0.3162}     {+0.0000 +1.0000 +0.0000}     {-0.3162 +0.0000 +0.9487}     {  O +Y  O };   # Fiber
#      14                 dispBeamColumn     {+0.9487 +0.0000 +0.3162}     {+0.0000 +1.0000 +0.0000}     {-0.3162 +0.0000 +0.9487}     {  O +Y  O };   # Fiber
