import numpy as np
from XCOM import mu_tot

### File Input Parameters

Z_range = np.array([1, 6, 13, 20, 26, 32, 40, 47, 55, 64, 74, 82, 92, 101, 102, 103, 104, 105, 106, 107])
E_range = ["4", "6", "10"]
lmbda_spacing = 10  # lambda spacing between simulations
N0 = 1e6            # thin target num_particles
N1 = 2146483647     # thick target num_particles
max_error = 5e-4
xml_path = "/nfs/home2/plalor/grasshopper/xml/gdml.xsd"

### Loading files to approximate the appropriate number of MC particles to run

path = "/Users/peter/Work/semiempirical_transparency/data/"
D = np.load(path + "D.npy")
D2 = np.load(path + "D2.npy")
E = np.load(path + "E.npy")
phi_4 = np.load(path + "phi_4MeV.npy")

def calcRelError(lmbda_arr, Z_arr, phi):
    phi0 = phi.copy()
    for lmbda, Z in zip(lmbda_arr, Z_arr):
        mu = mu_tot(E, Z)
        m0 = np.exp(-mu * lmbda)
        phi0 *= m0
    d0 = np.dot(D, phi0)
    sigma0 = np.sqrt(np.dot(D2, phi0))
    return sigma0 / d0

### Defining materials files

material_name = ["G4_H", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_N", "G4_O", "G4_F", "G4_Ne", "G4_Na", "G4_Mg", "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_Ar", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr", "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_Br", "G4_Kr", "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd", "G4_In", "G4_Sn", "G4_Sb", "G4_Te", "G4_I", "G4_Xe", "G4_Cs", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd", "G4_Pm", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho", "G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf", "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au", "G4_Hg", "G4_Tl", "G4_Pb", "G4_Bi", "G4_Po", "G4_At", "G4_Rn", "G4_Fr", "G4_Ra", "G4_Ac", "G4_Th", "G4_Pa", "G4_U", "G4_Np", "G4_Pu", "G4_Am", "G4_Cm", "G4_Bk", "G4_Cf"]
material_Z = np.arange(1, len(material_name)+1)
materials = {Z: material for (Z, material) in zip(material_Z, material_name)}

# add compound materials with unique Z identifier materials
compound_Z = {}
compound_w = {}

materials[101] = "G4_POLYETHYLENE"
compound_Z[101] = np.array([1, 6])
compound_w[101] = np.array([0.143711, 0.856289])

materials[102] = "G4_ALUMINUM_OXIDE"
compound_Z[102] = np.array([8, 13])
compound_w[102] = np.array([0.470749, 0.529251])

materials[103] = "G4_SILVER_CHLORIDE"
compound_Z[103] = np.array([17, 47])
compound_w[103] = np.array([0.247368, 0.752632])

materials[104] = "G4_LITHIUM_IODIDE"
compound_Z[104] = np.array([3, 53])
compound_w[104] = np.array([0.051858, 0.948142])

materials[105] = "G4_CADMIUM_TUNGSTATE"
compound_Z[105] = np.array([8, 48, 74])
compound_w[105] = np.array([0.177644, 0.312027, 0.510329])

materials[106] = "G4_GLASS_LEAD"
compound_Z[106] = np.array([8, 14, 22, 33, 82])
compound_w[106] = np.array([0.156453, 0.080866, 0.008092, 0.002651, 0.751938])

materials[107] = "G4_URANIUM_OXIDE"
compound_Z[107] = np.array([8, 92])
compound_w[107] = np.array([0.118502, 0.881498])

### Creating files

N_total = 0
SLURM_ARRAY_TASK_ID = 0
for E0 in E_range:
    phi = np.load(path + "phi_%sMeV.npy" % E0)
    for Z in Z_range:
        material = materials[Z]
        lmbda = 0
        while True:
            lmbda += lmbda_spacing
            if Z > 100:
                lmbda_arr = lmbda * compound_w[Z]
                Z_arr = compound_Z[Z]
            else:
                lmbda_arr = [lmbda]
                Z_arr = [Z]
            if calcRelError(lmbda_arr, Z_arr, phi_4)/np.sqrt(N1) > max_error:
                break
            if (Z == 1 or Z == 101) and 1.4*calcRelError(lmbda_arr, Z_arr, phi_4)/np.sqrt(N1) > max_error:
                break   ### hydrogen is slow
                    
            error = calcRelError(lmbda_arr, Z_arr, phi)
            N = int(N0 + (error / max_error)**2)
            assert N <= 2147483647
            
            filename = f"ID={SLURM_ARRAY_TASK_ID}-E={E0}MeV-lmbda={lmbda}-Z={Z}-N={N}.gdml"
            filestring = f"""<?xml version="1.0" encoding="UTF-8" standalone="no" ?>

<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="{xml_path}">
  
  <define>    
    <!-- target container -->
    <constant name="target_lambda" value="{lmbda}"/> <!-- target area density, g/cm^2 -->
    <constant name="target_x" value="50"/> <!-- target thickness in cm -->
    <constant name="target_y" value="50"/> <!-- target depth in cm -->
    <constant name="target_z" value="100"/> <!-- target height in cm -->
    <constant name="target_dist" value="100"/> <!-- distance from source to target in cm -->
    <constant name="target_rho" value="target_lambda/target_y"/> <!-- target density -->
    
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
    <quantity name="ProductionLowLimit" type="threshold" value="80" unit="keV" /> <!-- for neutron processes anything >1keV causes things to hang...set this to a high value for any other process to optimize computational time.  There are still some intricacies with this.  With high enough energy, rather than generating secondaries, all the energy loss will get tagged to the EnergyDeposited for the main particle.  So, the energy scoring (as determined by LighProducingParticle above) needs to be adjusted accordingly. -->
  </define>

<!-- OUTPUT FILTERS.  What data/entries do you want to keep? -->
  <define>
    <constant name="SaveSurfaceHitTrack" value="0"/> <!-- save entries which contain hit information, e.g. if you want to simulate the flux of particles -->
    <constant name="SaveTrackInfo" value="1"/> <!-- save individual track info (within an event).  This is useful for studying the physics of the various interactions -->
    <constant name="SaveEdepositedTotalEntry" value="0"/> <!--save entries which summarize the total deposited energy, e.g. in detector response simulations -->
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
    
    <!-- target -->
    <box name = "target_box" x="target_x" y="target_y" z="target_z" lunit= "cm"/>
    
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
    
    <!-- target -->
    <volume name="target_log">
      <materialref ref="Target"/>
      <solidref ref="target_box"/>
    </volume>
    
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
      
      <!-- target -->
      <physvol name="target_phys">
        <volumeref ref="target_log"/>
        <position name="target_pos" unit="cm" x="0" y="target_dist+target_y/2" z="target_z/2"/>
      </physvol>
            
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
