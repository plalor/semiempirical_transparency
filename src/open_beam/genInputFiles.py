import numpy as np
from XCOM import mu_tot

### File Input Parameters

E_range = ["4", "6", "10"]
max_error = 1e-5
xml_path = "/nfs/home2/plalor/grasshopper/xml/gdml.xsd"

### Loading files to approximate the appropriate number of MC particles to run

path = "/Users/peter/Work/semiempirical_transparency/data/"
D = np.load(path + "D.npy")
D2 = np.load(path + "D2.npy")
E = np.load(path + "E.npy")

def calcRelError(phi):
    d0 = np.dot(D, phi)
    sigma0 = np.sqrt(np.dot(D2, phi))
    return sigma0 / d0

### Creating files

N_total = 0
SLURM_ARRAY_TASK_ID = 0
for E0 in E_range:
    phi = np.load(path + "phi_%sMeV.npy" % E0)      
    error = calcRelError(phi)
    N = int((error / max_error)**2)
    assert N <= 2147483647 # Largest value for G4Int
        
    filename = f"ID={SLURM_ARRAY_TASK_ID}-E={E0}MeV-N={N}.gdml"
    filestring = f"""<?xml version="1.0" encoding="UTF-8" standalone="no" ?>

<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="{xml_path}">
  
  <define>    
    <!-- collimator properties -->
    <constant name="collimator_x" value="20"/>  <!-- length in cm -->
    <constant name="collimator_y" value="10"/>  <!-- width in cm -->
    <constant name="collimator_z" value="100"/> <!-- height in cm -->
    <constant name="collimator_sep" value="1"/> <!-- separation between collimators in cm -->
    
    <!-- detector properties -->
    <constant name="detector_x" value="1.0"/> <!-- length in cm -->
    <constant name="detector_y" value="3.0"/> <!-- width in cm -->
    <constant name="detector_z" value="1.0"/>  <!-- height in cm -->
    <constant name="detector_dist" value="200"/> <!-- distance from source to detector in cm -->
    
    <!-- do some calculations to determine these quantities -->
    <constant name="n_detectors" value="100"/> <!-- collimator_z // detector_z -->
    
    <!-- loop index -->
    <variable name="i"/>
  </define>
       
  <!-- THE OUTPUT -->
  <define>
    <constant name="TextOutputOn" value="1"/> <!-- the value should be either 1 (true) or 0 -->
    <constant name="BriefOutputOn" value="0"/> <!-- enable this if you want the shorter version of the text output -->
    <constant name="VRMLvisualizationOn" value="0"/> <!-- 1 means that you want a VRML file -->
    <constant name="EventsToAccumulate" value="50"/> <!-- number of tracks to accumulate in the visualization -->
  </define>
  
  <!-- CUTS...apply various cuts to make the computation more efficient -->
  <define>
    <constant name="LightProducingParticle" value="0"/> <!-- the particle which is actually producing light in the detector.  0 means ALL.  It will also kill all particles other than LightProducingParticle in the detector.  If in doubt set to 0. -->
    <constant name="LowEnergyCutoff" value="0."/><!-- The low energy cuttoff, MeV, for the main track. If in doubt set it to 0 -->
    <constant name="KeepOnlyMainParticle" value="0"/> <!-- if 1, the simulation will track only the main particle, as defined by ParticleNumber in the beam definition, OUTSIDE the detector volume.  For example, you'll need to set this to 0 to simulate bremmstrahlung, but to 1 for any transmission simulation. If in doubt, set to 0.-->
    <quantity name="ProductionLowLimit" type="threshold" value="80" unit="keV" /> <!-- for neutron processes anything >1keV causes things to hang...set this to a high value for any other process to optimize computational time.  There are still some intricacies with this.  With high enough energy, rather than generating secondaries, all the energy loss will get tagged to the EnergyDeposited for the main particle.  So, the energy scoring (as determined by LighProducingParticle above) needs to be adjusted accordingly. -->
  </define>

<!-- OUTPUT FILTERS.  What data/entries do you want to keep? -->
  <define>
    <constant name="SaveSurfaceHitTrack" value="0"/> <!-- save entries which contain hit information, e.g. if you want to simulate the flux of particles -->
    <constant name="SaveTrackInfo" value="0"/> <!-- save individual track info (within an event).  This is useful for studying the physics of the various interactions -->
    <constant name="SaveEdepositedTotalEntry" value="1"/> <!--save entries which summarize the total deposited energy, e.g. in detector response simulations -->
  </define>

  <!-- THE BEAM -->
  <define>
    <constant name="RandomGenSeed" value="1"/>
    <quantity name="BeamOffsetX"  type="coordinate" value="0" unit="cm"/>
    <quantity name="BeamOffsetY"  type="coordinate" value="0" unit="cm"/>
    <quantity name="BeamOffsetZ"  type="coordinate" value="0" unit="cm"/>
    <quantity name="BeamSize" type="coordinate" value="-1" unit="mm"/>
    <quantity name="BeamEnergy" type="energy" value="-2" unit="MeV"/> <!-- this is in MeV --> <!-- a negative number prompts reading input_spectrum.txt -->
    <quantity name="PhiMin" value="pi/2-atan(0.5*collimator_sep/detector_dist)"/>
    <quantity name="PhiMax" value="pi/2+atan(0.5*collimator_sep/detector_dist)"/>
    <quantity name="ThetaMin" value="pi/2"/>
    <quantity name="ThetaMax" value="pi/2-atan(collimator_z/detector_dist)"/>
    <constant name="EventsToRun" value="{N}"/>
    <constant name="ParticleNumber" value="22"/>
    <!-- e- is 11, gamma is 22, neutron is 2112, proton is 2212, alpha is 1000020040 -->
 
  </define>

  <!-- definition of solid geometries -->
  <solids>
    <!-- world volume -->
    <box name="world_solid" x="20" y="20" z="20" lunit="m"/>
    
    <!-- collimators -->
    <box name="collimator1" x="collimator_x" y="collimator_y" z="collimator_z" lunit="cm"/>
    <box name="collimator2" x="collimator_x" y="collimator_y" z="collimator_z" lunit="cm"/>
    
    <!-- cadmium tungstate detectors -->
    <loop for="i" from="1" to="n_detectors" step="1">
      <box name = "CdWO4_detector[i]" x="detector_x" y="detector_y" z="detector_z" lunit= "cm"/>
    </loop>
  </solids>


  <!-- PUTTING IT ALL TOGETHER -->
  <structure>
       
    <!-- collimators -->
    <volume name="collimator1_log">
      <materialref ref="G4_Pb"/>
      <solidref ref="collimator1"/>
    </volume>
    
    <volume name="collimator2_log">
      <materialref ref="G4_Pb"/>
      <solidref ref="collimator2"/>
    </volume>

    <!-- CdWO4 detector -->
    <loop for="i" from="1" to="n_detectors" step="1">
      <volume name="CdWO4_detector_log[i]">
        <materialref ref="G4_CADMIUM_TUNGSTATE"/>
        <solidref ref="CdWO4_detector[i]"/>
      </volume>
    </loop>
    
    <!-- top level world volume with all geometry elements -->
    <volume name="world_log">
      <materialref ref="G4_AIR"/>
      <solidref ref="world_solid"/>  <!-- world_solid This should NEVER be changed -->
            
      <!-- collimators -->
      <physvol name="collimator1_phys">
        <volumeref ref="collimator1_log"/>
        <position name="collimator1_pos" x="(collimator_x+collimator_sep)/2" y="detector_dist-collimator_y/2" z="collimator_z/2" unit="cm"/>
      </physvol>
      
      <physvol name="collimator2_phys">
        <volumeref ref="collimator2_log"/>
        <position name="collimator2_pos" x="-(collimator_x+collimator_sep)/2" y="detector_dist-collimator_y/2" z="collimator_z/2" unit="cm"/>
      </physvol>
      
      <!-- detector -->
      <loop for="i" from="1" to="n_detectors" step="1">
        <physvol name="det_phys[i]">
          <volumeref ref="CdWO4_detector_log[i]"/>
          <position name="CdWO4_detector_pos[i]" x="0" y="detector_dist+detector_x/2" z="(2*i-1)*detector_z/2" unit="cm"/>
        </physvol>
      </loop>

    </volume>

  </structure>

  <setup name="Default" version="1.0">
    <world ref="world_log"/>
  </setup>
</gdml>
"""
    with open(filename, "w") as f:
        f.write(filestring)
    SLURM_ARRAY_TASK_ID += 1
    N_total += N

print("Complete, with %.2e total particles and %d input files" % (N_total, SLURM_ARRAY_TASK_ID))
