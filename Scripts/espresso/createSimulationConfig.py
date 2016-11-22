#!/usr/bin/python

import os
from numpy import *

descr = "densePIL_eta0.3_elc-wall-charged"
subdescr = "extForceEqui01"
output_folder = "./generatedConfigs/"+descr+"_"+subdescr+"/"
tclScriptName = "polymer_sidechains_npt_wall_elc.tcl"
espresso_path = "~/espresso-master/"
pil_simulation_path = "/work/aweyman/pil_simulation/"
NNodes = 1 #no limit if 0
NCPUs = 16
maxComputationTime = "47:59:00"

#createSidechains = 1
N_P = 30
backboneFixed = 0
simulation_max_time = "150000" #use "+" for additional time units to simulation, without for absolute maximum simulation time
equi_dt = 0.01
equi_temp = 2.0 
enableCoulombInt = 1 #if 0, no p3m will be turned on during the simulation
coulombWarmup = 1
coulombWarmup_changeBjerrLenFactor = 0.1 #the smallest bjerrum length is this factor times l_bjerrum

vmd_writeCoordinates = 0
vmd_writeCoordinatesInterval = 1
vmd_writeHeader = 0

writeVelocities = 0
writeVelocitiesInterval = 1
velocity_writeHeader = 0

externalElField = 0.05 #0.05 #defines external force = externalElField*valency, ignored in tcl if 0

NPTEqui = 0
NPTPressure = 0.0036 #some artificial value in the order of 1 bar
NPTGamma0 = 1.0
NPTInitPistonMass = 1e-6
NPTMinPistonMass = 1e-6
NPTMaxPistonMass = 1
NPTPressureRelativeErrorTolerance = 0.05
NPTRecalculatePistonMassIterationInterval = 100
NPTRequiredSuccessfulIterations = 100

adjustVolume_exitAfterTargetValue = 0
adjustVolume_changeBoxlFactor = 0.0001
adjustVolume_p3mTuneIterationInterval = 100

checkpoint_iteration_interval = 200
save_all_checkpoints = 1 #if 1, every checkpoint for every checkpoint_iteration_interval (and for simulation abort) will be copied to a checkpointing subdir

calculateIntegratedCIDistribution = 0
distribution_iteration_interval = 1
distr_bins = 100
distr_log_flag = 1
distr_int_flag = 1
write_all_distr_values = 0

calculateFormFactor = 0
formfactor_iteration_interval = 10
formfac_bins = 100
formfac_q_min_fac = 0.01
formfac_q_max_fac = 100.0
write_all_formfac_values = 0

calculateRDF = 0
rdf_iteration_interval = 1
rdf_bins = 100
write_all_rdf_values = 0

calculateStructureFactor = 0
structurefactor_iteration_interval = 10
strufac_order = 5
write_all_strufac_values = 0

wallConstraint = 1
wall_LJ_eps = 1.75 #wall - counterions, WCA
wall_LJ_cut_factor = round(2**(1./6.),9) #wall - counterions, WCA
wall_LJ2_eps = 1.75 #wall - all types excl. CI
wall_LJ2_cut_factor = round(2**(1./6.),9) #2.5 #round(2**(1./6.),9) #wall - all types excl. CI
wall_charge_density = 0.05 #charges/sigma**2, typ. 0.05

enableELC = 1 #0 or 1
elcGapRatio = 0.2 #>0.0, typ. 0.2 (h/L ratio, with L+h=Lz)


#create folder if necessary
if not os.path.isdir(output_folder):
  os.makedirs(output_folder)

sigma = 1.0


#densePIL_eta0.3_v0x
MPC = [30]
l_B = [7.0]
#r = 0.0140625 #corresponds to box_l = 40 for N=30, MPC=30
charge_distance = [1]
sidechain_length = [0,1,2,3]


def calcRhoFromEta(deltaS, eta=0.3, decimals=8):
    sigma = 1.0
    return round(6.*eta/(pi*(deltaS+2)*sigma**3.), decimals)


submit_file = open(output_folder+"{0}_submit_{1}.sh".format(descr, subdescr), "w")

filecounter = 0

for sclen in sidechain_length:
  for N in MPC:
    for lb in l_B:    
      for dist in charge_distance:
        #write config file
        #conf_file = open(output_folder+"{0}_N{1}_rho{2}_lB{3}_dist{4}_simulation_config.txt".format(descr, N, r, lb, dist), "w")
        
        r = calcRhoFromEta(sclen) #chose accordingly for dense (calculation from rho) or dilute (fixed rho) systems
        
        conf_file = open(output_folder+"{0}_N{1}_rho{2}_lB{3}_dist{4}_sclen{5}_simulation_config.txt".format(descr, N, r, lb, dist, sclen), "w")
        
        conf_file.write("{tclvariable \n")
        conf_file.write("\t{{N_P {0}}}\n".format(N_P))
        conf_file.write("\t{{MPC {0}}}\n".format(N))
        #conf_file.write("\t{{createSidechains {0}}}\n".format(createSidechains))
        conf_file.write("\t{{createSidechains {0}}}\n".format(sclen))
        conf_file.write("\t{{backboneFixed {0}}}\n".format(backboneFixed))
        conf_file.write("\t{{calculateIntegratedCIDistribution {0}}}\n".format(calculateIntegratedCIDistribution))
        conf_file.write("\t{{distribution_iteration_interval {0}}}\n".format(distribution_iteration_interval))
        conf_file.write("\t{{distr_bins {0}}}\n".format(distr_bins))
        conf_file.write("\t{{distr_log_flag {0}}}\n".format(distr_log_flag))
        conf_file.write("\t{{distr_int_flag {0}}}\n".format(distr_int_flag))
        conf_file.write("\t{{write_all_distr_values {0}}}\n".format(write_all_distr_values))
        conf_file.write("\t{{calculateStructureFactor {0}}}\n".format(calculateStructureFactor))
        conf_file.write("\t{{structurefactor_iteration_interval {0}}}\n".format(structurefactor_iteration_interval))
        conf_file.write("\t{{strufac_order {0}}}\n".format(strufac_order))
        conf_file.write("\t{{write_all_strufac_values {0}}}\n".format(write_all_strufac_values))
        conf_file.write("\t{{calculateFormFactor {0}}}\n".format(calculateFormFactor))
        conf_file.write("\t{{formfactor_iteration_interval {0}}}\n".format(formfactor_iteration_interval))
        conf_file.write("\t{{formfac_bins {0}}}\n".format(formfac_bins))
        conf_file.write("\t{{formfac_q_min_fac {0}}}\n".format(formfac_q_min_fac))
        conf_file.write("\t{{formfac_q_max_fac {0}}}\n".format(formfac_q_max_fac))
        conf_file.write("\t{{write_all_formfac_values {0}}}\n".format(write_all_formfac_values))
        conf_file.write("\t{{calculateRDF {0}}}\n".format(calculateRDF))
        conf_file.write("\t{{rdf_iteration_interval {0}}}\n".format(rdf_iteration_interval))
        conf_file.write("\t{{rdf_bins {0}}}\n".format(rdf_bins))
        conf_file.write("\t{{write_all_rdf_values {0}}}\n".format(write_all_rdf_values))
        conf_file.write("\t{{monomer_density {0}}}\n".format(r))
        conf_file.write("\t{{side_bead_distance {0}}}\n".format(dist))
        conf_file.write("\t{{wallConstraint {0}}}\n".format(wallConstraint))
        conf_file.write("\t{{wall_LJ_eps {0}}}\n".format(wall_LJ_eps))
        conf_file.write("\t{{wall_LJ_cut_factor {0}}}\n".format(wall_LJ_cut_factor))
        conf_file.write("\t{{wall_LJ2_eps {0}}}\n".format(wall_LJ2_eps))
        conf_file.write("\t{{wall_LJ2_cut_factor {0}}}\n".format(wall_LJ2_cut_factor))
        conf_file.write("\t{{wall_charge_density {0}}}\n".format(wall_charge_density))
        conf_file.write("\t{{enableELC {0}}}\n".format(enableELC))
        conf_file.write("\t{{elcGapRatio {0}}}\n".format(elcGapRatio))
        conf_file.write("\t{LJ_eps 1.0}\n")
        conf_file.write("\t{LJ_cut_factor 1.122462048}\n")
        conf_file.write("\t{LJ2_eps 1.75}\n")
        conf_file.write("\t{LJ2_cut_factor 2.5}\n")
        conf_file.write("\t{fene_k 7.0}\n")
        conf_file.write("\t{fene_r 2.0}\n")
        conf_file.write("\t{{l_bjerrum {0}}}\n".format(lb))
        conf_file.write("\t{{enableCoulombInt {0}}}\n".format(enableCoulombInt))
        conf_file.write("\t{{coulombWarmup {0}}}\n".format(coulombWarmup))
        conf_file.write("\t{{coulombWarmup_changeBjerrLenFactor {0}}}\n".format(coulombWarmup_changeBjerrLenFactor))
        conf_file.write("\t{{externalElField {0}}}\n".format(externalElField))
        conf_file.write("\t{{equi_dt {0}}}\n".format(equi_dt))
        conf_file.write("\t{{equi_temp {0}}}\n".format(equi_temp))
        conf_file.write("\t{equi_gamma 1.0}\n")
        conf_file.write("\t{{vmd_writeCoordinates {0}}}\n".format(vmd_writeCoordinates))
        conf_file.write("\t{{vmd_writeHeader {0}}}\n".format(vmd_writeHeader))
        conf_file.write("\t{{vmd_writeCoordinatesInterval {0}}}\n".format(vmd_writeCoordinatesInterval))
        conf_file.write("\t{{writeVelocities {0}}}\n".format(writeVelocities))
        conf_file.write("\t{{velocity_writeHeader {0}}}\n".format(velocity_writeHeader))
        conf_file.write("\t{{writeVelocitiesInterval {0}}}\n".format(writeVelocitiesInterval))
        conf_file.write("\t{{checkpoint_iteration_interval {0}}}\n".format(checkpoint_iteration_interval))
        conf_file.write("\t{{save_all_checkpoints {0}}}\n".format(save_all_checkpoints))
        conf_file.write("\t{{NPTEqui {0}}}\n".format(NPTEqui))
        conf_file.write("\t{{NPTPressure {0}}}\n".format(NPTPressure))
        conf_file.write("\t{{NPTGamma0 {0}}}\n".format(NPTGamma0))
        conf_file.write("\t{{NPTInitPistonMass {0}}}\n".format(NPTInitPistonMass))
        conf_file.write("\t{{NPTMinPistonMass {0}}}\n".format(NPTMinPistonMass))
        conf_file.write("\t{{NPTMaxPistonMass {0}}}\n".format(NPTMaxPistonMass))
        conf_file.write("\t{{NPTPressureRelativeErrorTolerance {0}}}\n".format(NPTPressureRelativeErrorTolerance))
        conf_file.write("\t{{NPTRecalculatePistonMassIterationInterval {0}}}\n".format(NPTRecalculatePistonMassIterationInterval))
        conf_file.write("\t{{NPTRequiredSuccessfulIterations {0}}}\n".format(NPTRequiredSuccessfulIterations))
        conf_file.write("\t{{adjustVolume_exitAfterTargetValue {0}}}\n".format(adjustVolume_exitAfterTargetValue))
        conf_file.write("\t{{adjustVolume_changeBoxlFactor {0}}}\n".format(adjustVolume_changeBoxlFactor))
        conf_file.write("\t{{adjustVolume_p3mTuneIterationInterval {0}}}\n".format(adjustVolume_p3mTuneIterationInterval))
        conf_file.write("}")
        
        conf_file.close()
      
      
        #write slurm script file
        slurm_file = open(output_folder+"slurm_{0}_N{1}_rho{2}_lB{3}_dist{4}_sclen{5}.sh".format(descr, N, r, lb, dist, sclen), "w")
        
        slurm_file.write("#!/bin/bash\n")
        slurm_file.write("#SBATCH -J {0}_N{1}_rho{2}_lB{3}_dist{4}_sclen{5}\n".format(descr, N, r, lb, dist, sclen))
        if NNodes:
            slurm_file.write("#SBATCH -N {0}\n".format(NNodes))
        slurm_file.write("#SBATCH -n {0}\n".format(NCPUs))
        slurm_file.write("#SBATCH -x compute-big-0-0,compute-big-0-1\n")
	#slurm_file.write("#SBATCH --cpu_bind=socket")
	slurm_file.write("#SBATCH --time={0}\n".format(maxComputationTime))
        slurm_file.write("#SBATCH --mail-type=ALL\n")
        slurm_file.write("#SBATCH --mail-user=aweyman@icp.uni-stuttgart.de\n")
        #slurm_file.write("cd $SLURM_SUBMIT_DIR\n")
        #slurm_file.write("cd ~/pil_simulation\n")
        slurm_file.write("cd {0}\n".format(pil_simulation_path))
        slurm_file.write("mpirun {0}Espresso {8} {1}_N{2}_rho{3}_lB{4}_dist{5}_sclen{6} {7}".format(espresso_path, descr, N, r, lb, dist, sclen, simulation_max_time, tclScriptName))
	#slurm_file.write("{0}Espresso polymer_sidechains.tcl {1}_N{2}_rho{3}_lB{4}_dist{5}_sclen{6} {7}".format(espresso_path, descr, N, r, lb, dist, sclen, simulation_max_time))
      
        slurm_file.close()
      
        filecounter+=1
      
      
        #write slurm submit file
        submit_file.write("sbatch slurm_{0}_N{1}_rho{2}_lB{3}_dist{4}_sclen{5}.sh\n".format(descr, N, r, lb, dist, sclen))
        submit_file.write("sleep 1\n")


submit_file.close()

print "{0} configurations created.".format(filecounter)
print "Wrote all data to {0}".format(output_folder)