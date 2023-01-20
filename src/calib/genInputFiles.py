import numpy as np
from XCOM import mu_tot

### File Input Parameters

zRange = np.array([6, 26, 82])
lmbdaRange = np.array([50, 100, 150])
num_jobs = 20
max_error = 4e-5
#xml_path = "/Users/peter/Work/grasshopperPeter/xml/gdml.xsd"
xml_path = "/home/plalor/grasshopperPeter/xml/gdml.xsd"

### Loading files to approximate the appropriate number of MC particles to run

path = "/Users/peter/Work/radiography/data/"
R = np.load(path + "R_10.npy")
E_g = np.load(path + "E_g_10.npy")
E_dep = np.load(path + "E_dep_10.npy")
b_4 = np.load(path + "b4MeV_10.npy")

def calcRelError(lmbda, Z, b):
    atten = mu_tot(E_g, Z)
    m0 = np.exp(-atten * lmbda)
    d0 = np.dot(R.T @ E_dep, m0 * b)
    sigma0 = np.sqrt(np.dot(R.T @ E_dep**2, m0 * b))
    return sigma0 / d0

### Defining materials files

material_name = ["G4_H", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_N", "G4_O", "G4_F", "G4_Ne", "G4_Na", "G4_Mg", "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_Ar", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr", "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_Br", "G4_Kr", "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd", "G4_In", "G4_Sn", "G4_Sb", "G4_Te", "G4_I", "G4_Xe", "G4_Cs", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd", "G4_Pm", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho", "G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf", "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au", "G4_Hg", "G4_Tl", "G4_Pb", "G4_Bi", "G4_Po", "G4_At", "G4_Rn", "G4_Fr", "G4_Ra", "G4_Ac", "G4_Th", "G4_Pa", "G4_U", "G4_Np", "G4_Pu", "G4_Am", "G4_Cm", "G4_Bk", "G4_Cf"]
material_Z = np.arange(1, len(material_name)+1)
materials = {Z: material for (Z, material) in zip(material_Z, material_name)}

### Creating files

for E in ["10.3", "5.5"]:
    b = np.load(path + "b%sMeV_10.npy" % E)
    for Z in zRange:
        material = materials[Z]
        for lmbda in lmbdaRange:
            error = calcRelError(lmbda, Z, b)
            N = int((error / max_error)**2 / num_jobs)
            
            filename = "E=%sMeV-lmbda=%d-Z=%d-N=%d.gdml" % (E, lmbda, Z, N)
            filestring = f"""<?xml version="1.0" encoding="UTF-8" standalone="no" ?>

<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="{xml_path}">
  
<define>
  <!-- target material -->
  <constant name="target_lambda" value="{lmbda}"/> <!-- target area density, g/cm^2 -->
  <constant name="x_target" value="100"/> <!-- target thickness in cm -->
  <constant name="y_target" value="200"/> <!-- target depth in cm -->
  <constant name="z_target" value="400"/> <!-- target height in cm -->
  <constant name="dist_target" value="350"/> <!-- distance from source to target in cm -->
  <constant name="target_rho" value="target_lambda/x_target"/> <!-- target density -->
  <!-- collimator properties -->
  <constant name="x_collimator" value="10"/>  <!-- thickness in cm -->
  <constant name="y_collimator" value="20"/>  <!-- depth in cm -->
  <constant name="z_collimator" value="400"/> <!-- height in cm -->
  <constant name="sep_collimator" value="1"/> <!-- separation between collimators in cm -->
  <!-- detector properties -->
  <constant name="x_detector" value="3.0"/> <!-- thickness in cm -->
  <constant name="y_detector" value="1.5"/> <!-- depth in cm -->
  <constant name="z_detector" value="0.46"/>  <!-- height in cm -->
  <constant name="dist_detector" value="700"/> <!-- distance from source to detector in cm -->
  <!-- do some calculations to determine these quantities -->
  <constant name="n_detectors" value="869"/> <!-- z_collimator // z_detector = 869 -->
  <!-- loop index (for generating detector stack -->
  <variable name="i"/>
</define>
   
   <materials>
    
    <material name="Target" state="solid">
      <D unit="g/cm3" value="target_rho"/>
      <fraction n="1" ref="{material}"/>
    </material>
    
  </materials>
    
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
    <quantity name="ProductionLowLimit" type="threshold" value="100" unit="keV" /> <!-- for neutron processes anything >1keV causes things to hang...set this to a high value for any other process to optimize computational time.  There are still some intricacies with this.  With high enough energy, rather than generating secondaries, all the energy loss will get tagged to the EnergyDeposited for the main particle.  So, the energy scoring (as determined by LighProducingParticle above) needs to be adjusted accordingly. -->
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
    <quantity name="BeamEnergy" type="energy" value="-1" unit="MeV"/> <!-- this is in MeV --> <!-- a negative number prompts reading input_spectrum.txt -->
    <quantity name="PhiMin" value="-atan(0.5*sep_collimator/dist_detector)"/>
    <quantity name="PhiMax" value="atan(0.5*sep_collimator/dist_detector)"/>
    <quantity name="ThetaMin" value="atan(dist_detector/z_collimator)"/>
    <quantity name="ThetaMax" value="pi/2"/>
    <constant name="EventsToRun" value="{N}"/>
    <constant name="ParticleNumber" value="22"/>
    <!-- e- is 11, gamma is 22, neutron is 2112, proton is 2212, alpha is 1000020040 -->
 
  </define>

  <!-- definition of solid geometries -->
  <solids>
    <!-- world volume -->
    <box name="world_solid" x="20" y="20" z="20" lunit="m"/>
    
    <!-- target -->
    <box name = "target_box" x="x_target" y="y_target" z="z_target" lunit= "cm"/>
    
    <!-- collimators -->
    <box name = "collimator_1" x="x_collimator" y="y_collimator" z="z_collimator" lunit= "cm"/>
    <box name = "collimator_2" x="x_collimator" y="y_collimator" z="z_collimator" lunit= "cm"/>
    
    <!-- cadmium tungstate detectors -->
    <loop for="i" from="1" to="n_detectors" step="1">
      <box name = "CdWO4_detector[i]" x="x_detector" y="y_detector" z="z_detector" lunit= "cm"/>
    </loop>
  </solids>


  <!-- PUTTING IT ALL TOGETHER -->
  <structure>
       
     <!-- target -->
     <volume name="target_log">
       <materialref ref="Target"/>
       <solidref ref="target_box"/>
     </volume>
       
     <!-- collimators -->
     <volume name="collimator_1_log">
       <materialref ref="G4_Pb"/>
       <solidref ref="collimator_1"/>
     </volume>
     
     <volume name="collimator_2_log">
       <materialref ref="G4_Pb"/>
       <solidref ref="collimator_2"/>
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
      
      <!-- target -->
      <physvol name="target_phys">
        <volumeref ref="target_log"/>
        <position name="target_pos" unit="cm" x="dist_detector/2" y="0" z="z_target/2"/>
      </physvol>
      
      <!-- collimators -->
      <physvol name="collimator_1_phys">
        <volumeref ref="collimator_1_log"/>
        <position name="collimator_1_pos" unit="cm" x="dist_detector-x_collimator/2" y="(y_collimator+sep_collimator)/2" z="z_collimator/2"/>
      </physvol>
      
      <physvol name="collimator_2_phys">
        <volumeref ref="collimator_2_log"/>
        <position name="collimator_2_pos" unit="cm" x="dist_detector-x_collimator/2" y="-(y_collimator+sep_collimator)/2" z="z_collimator/2"/>
      </physvol>
      
      <!-- detector -->
      <loop for="i" from="1" to="n_detectors" step="1">
        <physvol name="det_phys[i]">
          <volumeref ref="CdWO4_detector_log[i]"/>
          <position name="CdWO4_detector_pos[i]" unit="cm" x="dist_detector+x_detector/2" y="0" z="(2*i-1)*z_detector/2"/>
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

    ### Generating file with no target
      
    error = calcRelError(0, 26, b)
    N = int((error / max_error)**2)
        
    filename = "E=%sMeV-lmbda=0-N=%d.gdml" % (E, N)
    filestring = f"""<?xml version="1.0" encoding="UTF-8" standalone="no" ?>

<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="{xml_path}">

<define>
  <!-- collimator properties -->
  <constant name="x_collimator" value="10"/>  <!-- thickness in cm -->
  <constant name="y_collimator" value="20"/>  <!-- depth in cm -->
  <constant name="z_collimator" value="400"/> <!-- height in cm -->
  <constant name="sep_collimator" value="1"/> <!-- separation between collimators in cm -->
  <!-- detector properties -->
  <constant name="x_detector" value="3.0"/> <!-- thickness in cm -->
  <constant name="y_detector" value="1.5"/> <!-- depth in cm -->
  <constant name="z_detector" value="0.46"/>  <!-- height in cm -->
  <constant name="dist_detector" value="700"/> <!-- distance from source to detector in cm -->
  <!-- do some calculations to determine these quantities -->
  <constant name="n_detectors" value="869"/> <!-- z_collimator // z_detector = 869 -->
  <!-- loop index (for generating detector stack -->
  <variable name="i"/>
</define>

   <materials>

  </materials>

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
    <quantity name="ProductionLowLimit" type="threshold" value="100" unit="keV" /> <!-- for neutron processes anything >1keV causes things to hang...set this to a high value for any other process to optimize computational time.  There are still some intricacies with this.  With high enough energy, rather than generating secondaries, all the energy loss will get tagged to the EnergyDeposited for the main particle.  So, the energy scoring (as determined by LighProducingParticle above) needs to be adjusted accordingly. -->
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
    <quantity name="BeamEnergy" type="energy" value="-1" unit="MeV"/> <!-- this is in MeV --> <!-- a negative number prompts reading input_spectrum.txt -->
    <quantity name="PhiMin" value="-atan(0.5*sep_collimator/dist_detector)"/>
    <quantity name="PhiMax" value="atan(0.5*sep_collimator/dist_detector)"/>
    <quantity name="ThetaMin" value="atan(dist_detector/z_collimator)"/>
    <quantity name="ThetaMax" value="pi/2"/>
    <constant name="EventsToRun" value="{N}"/>
    <constant name="ParticleNumber" value="22"/>
    <!-- e- is 11, gamma is 22, neutron is 2112, proton is 2212, alpha is 1000020040 -->

  </define>

  <!-- definition of solid geometries -->
  <solids>
    <!-- world volume -->
    <box name="world_solid" x="20" y="20" z="20" lunit="m"/>

    <!-- collimators -->
    <box name = "collimator_1" x="x_collimator" y="y_collimator" z="z_collimator" lunit= "cm"/>
    <box name = "collimator_2" x="x_collimator" y="y_collimator" z="z_collimator" lunit= "cm"/>

    <!-- cadmium tungstate detectors -->
    <loop for="i" from="1" to="n_detectors" step="1">
      <box name = "CdWO4_detector[i]" x="x_detector" y="y_detector" z="z_detector" lunit= "cm"/>
    </loop>
  </solids>


  <!-- PUTTING IT ALL TOGETHER -->
  <structure>

     <!-- collimators -->
     <volume name="collimator_1_log">
       <materialref ref="G4_Pb"/>
       <solidref ref="collimator_1"/>
     </volume>

     <volume name="collimator_2_log">
       <materialref ref="G4_Pb"/>
       <solidref ref="collimator_2"/>
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
      <physvol name="collimator_1_phys">
        <volumeref ref="collimator_1_log"/>
        <position name="collimator_1_pos" unit="cm" x="dist_detector-x_collimator/2" y="(y_collimator+sep_collimator)/2" z="z_collimator/2"/>
      </physvol>

      <physvol name="collimator_2_phys">
        <volumeref ref="collimator_2_log"/>
        <position name="collimator_2_pos" unit="cm" x="dist_detector-x_collimator/2" y="-(y_collimator+sep_collimator)/2" z="z_collimator/2"/>
      </physvol>

      <!-- detector -->
      <loop for="i" from="1" to="n_detectors" step="1">
        <physvol name="det_phys[i]">
          <volumeref ref="CdWO4_detector_log[i]"/>
          <position name="CdWO4_detector_pos[i]" unit="cm" x="dist_detector+x_detector/2" y="0" z="(2*i-1)*z_detector/2"/>
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
