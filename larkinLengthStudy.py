# The purpose of this file is to run specific modifications to noise profiles and simulation cells in MoDELib to study the length effects experienced by dislocations in multicomponent systems
# Authors: Hyunsoo Lee & Liam Myhill

import pandas as pd
import os

#define a function to clean the contents of a uniformLoadController
#return:nothing
def cleanUFL(ufl):

    currentDir=os.getcwd() #store the location of the currentDir
    
    os.chdir(ufl) #go to the ufl
    
    os.system('rm evl/* F/*') #clean the evl and F directories
    
    os.chdir(currentDir) #return to the original directory
    
    return()
    
#function to run the microstructureGenerator tool
#return:nothing
def MG(toolsDir,ufl):

    print('Running microstructureGenerator...')

    currentDir=os.getcwd()

    os.chdir(toolsDir+'/MicrostructureGenerator/build')
    os.system('./microstructureGenerator '+ufl)

    os.chdir(currentDir)
    return()

#function to run DDomp
#return:nothing
def DDomp(toolsDir,ufl):

    print('Running DDomp')

    currentDir=os.getcwd()

    os.chdir(toolsDir+'/DDomp/build')
    os.system('./DDomp '+ufl)

    os.chdir(currentDir)

    return()

#function to run DD2Ovito
#return:nothing
def DD2OvitoVTK(toolsDir,ufl):

    print('Running DD2OvitoVtk')

    currentDir=os.getcwd()

    os.chdir(toolsDir+'/DD2OvitoVtk/build')
    os.system('./DD2OvitoVtk '+ufl)

    os.chdir(currentDir)

    return()
    
#define a function to copy all configuration information directly to a storage location (add functionality to capture material file and microstructure file)
#return: nothing
def copyConfigAndResults(ufl,dataLoc):
    
    cwd=os.getcwd()
    
    os.chdir(ufl)
    
    os.system('cp -r inputFiles/ evl/ F/ '+str(dataLoc))
    
    os.chdir(cwd)
    
    return()
    
#function to adjust the constant stress applied in a MoDELib2 DDD sim
#return:nothing
def adjustConstStress(ufl,component,value):

    cwd=os.getcwd()

    inputFilesDir=ufl+'/inputFiles/'

    if component=='11':

        stressMat=[[value,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]

    elif component=='22':

        stressMat=[[0.0,0.0,0.0],[0.0,value,0.0],[0.0,0.0,0.0]]

    elif component=='33':

        stressMat=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,value]]

    elif component=='12' or component=='21':

        stressMat=[[0.0,value,0.0],[value,0.0,0.0],[0.0,0.0,0.0]]

    elif component=='13' or component=='31':

        stressMat=[[0.0,0.0,value],[0.0,0.0,0.0],[value,0.0,0.0]]

    elif component=='23' or component=='32':

        stressMat=[[0.0,0.0,0.0],[0.0,0.0,value],[0.0,value,0.0]]

    os.chdir(inputFilesDir)

    with open('uniformExternalLoadController.txt','r') as f:
        lines=f.readlines()
        
    string='ExternalStress0'

    for i in range(np.size(lines,axis=0)):
        string_loc_itr=re.search(string,lines[i])
        if string_loc_itr:
            string_loc=i
            
    lines[string_loc]=string+' ='+str(stressMat[0]).replace('[','').replace(']','').replace(',',' ')+'\n'+str(stressMat[1]).replace('[','').replace(']','').replace(',',' ')+'\n'+str(stressMat[2]).replace('[','').replace(']','').replace(',',' ')+';\n'

    with open('uniformExternalLoadController.txt','w') as f:
        for line in lines:
            f.write(line)

    os.chdir(cwd)

    return()
    
#define a function to determine whether or not the CRSS has been reached for a single run of MoDEL
#return(df)
def findCRSS(df,numSteps):
    regressionTol=3

    slopeTol=2e-4 #increase if CRSS plot doesn't demonstrate CRSS (normally 1e-4)

    df['slopeBP_12']=df['betaP_12'].rolling(regressionTol).apply(lambda s: linregress(s.reset_index())[0])

    df['slopeBP_12']=df['slopeBP_12'].fillna(0)

    if np.nanmean(df.loc[regressionTol+1:,'slopeBP_12'])>slopeTol and np.nanmean(df.loc[(4/5)*numSteps:,'slopeBP_12']>slopeTol):
        CRSS_reached=1
    else:
        CRSS_reached=0

    print('Total Avg: '+str(np.nanmean(df.loc[regressionTol+1:,'slopeBP_12'])))
    print('Final Avg: '+str(np.nanmean(df.loc[(4/5)*numSteps:,'slopeBP_12'])))
    print(df['slopeBP_12'])
    return(df,CRSS_reached)

#define a function to modify and run generatepolycrystal.m in inputFiles (currently tailored to needs of larkin length study) pass in diag as a string of 3 int values
#update to include rotations as necessary
#return:Nothing
def changeBoxSettings(ufl,diag):

    cwd=os.getcwd()

    os.chdir(ufl+'/inputFiles')
    
    with open('GeneratePolyCrystalFile.m','r') as f:
        lines=f.readlines()
        
    print(lines)
        
    string='Fs='
    regEx='^Fs=diag\(\[\s*\d+\s*\d+\s*\d+\s*\]\);\n$'
    
    for i in range(np.size(lines,axis=0)):
        FLoc_itr=re.search(regEx,lines[i])
        if FLoc_itr:
            stringLoc=i
            
            
    print(lines[stringLoc])
    
    lines[stringLoc]='Fs=diag(['+diag+']);'
    
    with open('GeneratePolyCrystalFile.m','w') as f:
        for line in lines:
            f.write(line)
            
    ### run the MATLAB script
    
    octave.addpath(ufl+'/inputFiles')
    
    octave.run('GeneratePolyCrystalFile.m')
        
    return()
    
#define a function to modify the noise profile calculated for MoDELib solid solution nosie
#return:noting
def changeNoise(seedNo,ufl,materialFile,MSSS_SI):

    print('Modifying the noise of MoDELib configuration...')

    currentDir=os.getcwd()

    os.chdir(ufl+'/inputFiles')

    string='seed='

    with open('polycrystal.txt','r') as polycrystal:
        lines=polycrystal.readlines()

    for i in range(np.size(lines,axis=0)):
        if string in lines[i]:
            stringLoc=i

    lines[stringLoc]='seed='+str(seedNo)+'; #stochastic seed No. for solid solution noise\n'

    with open('polycrystal.txt','w') as polycrystal:
        for i in range(np.size(lines,axis=0)):
            polycrystal.write(lines[i])
            
    DisDyn=ufl.rsplit('/',2)[0]+'/' #path to dislocationDynamics
    
    matLib=DisDyn+'MaterialsLibrary' #path to matlib

    os.chdir(matLib)

    string='MSSS_SI='

    with open(materialFile,'r') as f:
        lines=f.readlines()

    for i in range(np.size(lines,axis=0)):
        MSSSLoc_itr=re.search(string,lines[i])
        if MSSSLoc_itr:
            stringLoc=i

    lines[stringLoc]='MSSS_SI='+str(MSSS_SI)+'; #[Pa^2] Mean Square Shear Stress\n'

    with open(materialFile,'w') as f:
        for i in range(np.size(lines,axis=0)):
            f.write(lines[i])

    os.chdir(currentDir)

    return()

#define a function to extract material information
#return(mu,b,material,rho)
def getMaterialInfo(ufl,matLib):

    currentDir=os.getcwd()

    string='materialFile'

    os.chdir(ufl+'/inputFiles')
    with open('polycrystal.txt','r') as f:
        lines=f.readlines()

    for i in range(np.size(lines,axis=0)):
        material_Loc_itr=re.search(string,lines[i])
        if material_Loc_itr:
            matLine=lines[i]


    material=os.path.basename(matLine.split('/')[-1])   # returns Mat.txt
    material=material.strip('; \n')

    os.chdir(matLib)

    string1='b_SI'
    string2='mu0_SI'
    string3='rho_SI'

    with open(material,'r') as f:
        lines=f.readlines()

    for i in range(np.size(lines,axis=0)):
        bLoc_itr=re.search(string1,lines[i])
        muLoc_itr=re.search(string2,lines[i])
        rhoLoc_itr=re.search(string3,lines[i])
        if bLoc_itr:
            bLine=lines[i]
        if muLoc_itr:
            muLine=lines[i]
        if rhoLoc_itr:
            rhoLine=lines[i]


            
    b=float(bLine.strip('b_SI=').strip('; # [m] \t\tBurgers vector magnitude \n'))
    mu=float(muLine.strip('mu0_SI=').strip(';\t# [Pa] \t\ttemperature-independent shear modulus coeff in mu=mu0+mu1*T\n'))
    rho=float(rhoLine.strip('rho_SI=').strip(';\t# [kg/m^3]\tmass density\n'))
    return(mu,b,material,rho)
    
#define a function that runs a full convergence study on the CRSS obtained from MoDELib2
#return: nothing
def runCRSS_convergence(directory,startStress,endStress,origin,ssID,numNodesMat,exitFaceID,initialDipoleGlideSteps,dipoleHeight,partialsActive,numSimSteps,PBCs,ewaldLengthFactor,temperature,langevinThermostat,langevinSeed,quadPerLengthMat,numRelaxSteps,writeLoc):


    DisDyn=directory.rsplit('/',2)[0]+'/' #path to dislocationDynamics
    MoDELRoot=directory.rsplit('/',4)[0]+'/' #path to root
    toolsDir=MoDELRoot+'tools' #path to tools
    matLib=DisDyn+'MaterialsLibrary' #path to matlib
    micLib=DisDyn+'MicrostructureLibrary' #path to miclib
    DDsegLib=toolsDir+'/DDsegments/' #path to DDseg

    CRSS_dataMat=[[]]

    mu,b,material,rho=getMaterialInfo(directory,matLib)
    
    Fmat,V=getboxDim(directory)
    
    numTrials=50 #number of stress cases to run in order to detect CRSS
    
    stressRange_Pa=np.linspace(startStress,endStress,num=numTrials) #define the stress range
    
    stressRange_Mu=stressRange_Pa/mu #convert the stressRange to code units
    
    os.chdir(directory) #move to ufl
    
    os.system('rm evl/* F/*') #clean the evl&F directories
    
    runRelaxation(directory,numRelaxSteps) #specify the number of relaxation steps
                
    #    for stress in stressRange_Mu: #loop over stress values to determine CRSS
    
    for quadPerLength in quadPerLengthMat: #loop over prescribed quadPerLength values
    
        adjustSimulationSettings(directory,numSimSteps,PBCs,ewaldLengthFactor,temperature,langevinThermostat,langevinSeed,quadPerLength) #adjust settings for individual simulations
        
        for numNodes in numNodesMat: #loop over prescribed nodal values

            adjustDipoleMicrostructure(directory,material,origin,ssID,numNodes,exitFaceID,initialDipoleGlideSteps,dipoleHeight,partialsActive) #adjust settings in the microstructureFile
            
            count=0 #counts the number of failed CRSS detections
            
            CRSS_reached=0 # a value that changes to 1 when the CRSS has been reached
            
            while CRSS_reached==0:
            
                currentStress_Mu=stressRange_Mu[count]
                
                currentStress_Pa=currentStress_Mu*mu
                
                os.chdir(directory) #move to ufl
    
                os.system('rm evl/* F/*') #clean the evl&F directories
                
                adjustConstStress(directory,'12',currentStress_Mu) #change the stress
                    
                MG(toolsDir,directory) #run the microstructureGenerator to create the dipole
                    
                DDomp(toolsDir,directory) #run dynamics to evolve the trajectory
                        
                F_frame=pullF(directory,mu,rho,V,b,Fmat[2],partialsActive)

                F_frame, CRSS_reached = findCRSS(F_frame,numSimSteps)

                count+=1
                
                if count > np.size(stressRange_Mu,axis=0)-1:
                    print('The CRSS is greater than the maximum of the prescribed range')
                    break
                
            if CRSS_reached==0:
                #                print('The CRSS is greater than the maximum of the prescribed range')
                CRSS_dataMat.append([numNodes,quadPerLength,CRSS_reached,currentStress_Pa])
            else:
                print('The CRSS is detected as '+str(currentStress_Pa/1e6)+'MPa with '+str(numNodes)+' nodes and quadPerLength='+str(quadPerLength)+'\n')
                
                DD2OvitoVTK(toolsDir,ufl)
                
                
                
                CRSS_dataMat.append([numNodes,quadPerLength,CRSS_reached,currentStress_Pa])
    
        #    CRSS_dataMat=filter([],CRSS_dataMat)

    
    
    print('CRSS_DATA: \n')
    print(str(CRSS_dataMat))
    
    
    nodesData=[]
    quadPointData=[]
    CRSSData=[]
    
    #print CRSS_dataMat to the writeLoc
    
    os.chdir(writeLoc)
    
    with open('CRSSdata.txt','w') as f:
        for line in CRSS_dataMat:
            f.write(str(line))
    
    for data in CRSS_dataMat:
        if len(data)==0:
            continue
        if data[2]==0.0:
            continue
        else:
            nodesData.append(data[0])
            quadPointData.append(data[1])
            CRSSData.append(data[3])
        
    plot_with_colorbar(nodesData,CRSSData,quadPointData)
    
        #    mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        #
        #    # Using contourf to provide my colorbar info, then clearing the figure
        #    Z = [[0,0],[0,0]]
        #    step=1
        #    min, max=(0,1)
        #    levels = range(min,max+step,step)
        #    CS3 = plt.contourf(Z, levels, cmap=mymap)
        #    plt.clf()
        #
        #    X=nodesData
        #    Y=CRSSData
        #    Z=quadPointData
        #    for x,y,z in zip(X,Y,Z):
        #        # setting rgb color based on z normalized to my range
        #        r = (float(z)-min)/(max-min)
        #        g = 0
        #        b = 1-r
        #        plt.plot(x,y,color=(r,g,b))
        #    plt.colorbar(CS3) # using the colorbar info I got from contourf
        #    plt.savefig('CRSS_convergence.png')
        #    plt.show()

    return()

#function to extract box size information
#return([xdim,ydim,zdim],V)
def getboxDim(directory):
    os.chdir(directory+'/inputFiles')
    string='F='
    with open('polycrystal.txt') as f:
        lines=f.readlines()

    for i in range(np.size(lines,axis=0)):
        A_loc_itr=re.search(string,lines[i])
        if A_loc_itr:
            A_loc=i
    Amat=lines[A_loc:A_loc+3]
    for i in range(np.size(Amat,axis=0)):
            Amat[i]=Amat[i].split()
    xdim=float(re.findall("\d+\.\d+",Amat[0][0])[0])
    ydim=float(Amat[1][1])
    zdim=float(Amat[2][2])
    
    V=xdim*ydim*zdim
    # filter_object = list(filter(lambda a: string in a, lines)) #POTENTIALLY USEFUL LATER ON
    # filter_object_str=str(filter_object)
    # A_loc=lines.index(filter_object_str.replace('[','').replace(']',''))
    return([xdim,ydim,zdim],V)

#function to run relaxation steps
#return:nothing
def runRelaxation(ufl,numRelaxSteps):

    cwd=os.getcwd()
    
    inputFilesDir=ufl+'/inputFiles'
    
    os.chdir(inputFilesDir)
    
    with open('uniformExternalLoadController.txt','r') as f:
        lines=f.readlines()
        
    string='relaxSteps'
        
    for i in range(np.size(lines,axis=0)):
        string_loc_itr=re.search(string,lines[i])
        if string_loc_itr:
            string_loc=i
            
    lines[string_loc]=string+' ='+str(numRelaxSteps).replace('[','').replace(']','').replace(',',' ')+';\n'

    return()
    
#function to adjust all pertinent simulation settings
#return:nothing
def adjustSimulationSettings(ufl,numSimSteps,PBCs,ewaldLengthFactor,temperature,langevinThermostat,langevinSeed,quadPerLength):

    cwd=os.getcwd()
    
    inputFilesDir=ufl+'/inputFiles/'

    #Lines of DD.txt
    string1='periodicImageSize'
    string2='Nsteps'
    string3='EwaldLengthFactor'
    string4='use_stochasticForce'
    string5='stochasticForceSeed'
    
    #Lines of polycrystal.txt
    string6='absoluteTemperature'
    
    os.chdir(inputFilesDir)
    
    with open('DD.txt','r') as f:
        lines=f.readlines()
        
        for i in range(np.size(lines,axis=0)):
            string1_loc_itr=re.search(string1,lines[i])
            string2_loc_itr=re.search(string2,lines[i])
            string3_loc_itr=re.search(string3,lines[i])
            string4_loc_itr=re.search(string4,lines[i])
            string5_loc_itr=re.search(string5,lines[i])
            if string1_loc_itr:
                string1_loc=i
            elif string2_loc_itr:
                string2_loc=i
            elif string3_loc_itr:
                string3_loc=i
            elif string4_loc_itr:
                string4_loc=i
            elif string5_loc_itr:
                string5_loc=i
                
    lines[string1_loc]=string1+' ='+str(PBCs).replace('[','').replace(']','').replace(',',' ')+';\n'
    lines[string2_loc]=string2+' ='+str(numSimSteps).replace('[','').replace(']','').replace(',',' ')+';\n'
    lines[string3_loc]=string3+' ='+str(ewaldLengthFactor).replace('[','').replace(']','').replace(',',' ')+';\n'
    lines[string4_loc]=string4+' ='+str(langevinThermostat).replace('[','').replace(']','').replace(',',' ')+';\n'
    lines[string5_loc]=string5+' ='+str(langevinSeed).replace('[','').replace(']','').replace(',',' ')+';\n'
    
    with open('DD.txt','w') as f:
        for line in lines:
            f.write(line)
        
    with open('polycrystal.txt','r') as f:
        lines=f.readlines()
        
    for i in range(np.size(lines,axis=0)):
        string6_loc_itr=re.search(string6,lines[i])
        if string6_loc_itr:
            string6_loc=i
            
    lines[string6_loc]=string6+' ='+str(temperature).replace('[','').replace(']','').replace(',',' ')+';\n'
    
    with open('polycrystal.txt','w') as f:
        for line in lines:
            f.write(line)
    
    os.chdir(cwd)
    
    return()
    
# define a function to plot data with 3 entries
#return: nothing
def plot_with_colorbar(x, y, z):
    fig, ax = plt.subplots()
    scatter = ax.scatter(x, y, c=z, cmap='viridis')
    cbar = plt.colorbar(scatter)
    cbar.set_label('Quadrature Point Density')
    ax.set_xlabel('Nodal Density')
    ax.set_ylabel('CRSS')
    ax.set_title('Convergence Study')
    plt.savefig('CRSS_convergence.png')
    plt.show()
    ## Example usage:
    #x = np.random.rand(100)
    #y = np.random.rand(100)
    #z = np.random.rand(100) * 100  # Third vector with values for color
    #plot_with_colorbar(x, y, z)
    return()
    
def main():

    
    #point to necessary directories
    ufl=input('Point to the desired uniformLoadController: ')
    dataLoc=input('Point to the desired data storage location: ')

    #specify simulation settings
    seedMat=[1,2,3]
    crit_dim='x'
    boxLengthMat=[1000,700,500,300,100,50,25,10]
    MSSS_SI_Mat=[]
    
    #specify the shape of the simulation cell
    cubicCell=1
    tailoredCell=0
    
    

    #modify the noise profile on the glide plane
#    changeNoise(seedNo,ufl,materialFile,MSSS_SI)

    #modify the size of the dislocation cell
    for boxLength in boxLengthMat:
    
        cleanUFL(ufl)
        
        if tailoredCell==1 and cubicCell==0:
            diag=str(boxLength)+' 120 60'
        elif cubicCell==1 and tailoredCell==0:
            diag=str(boxLength)+' '+str(boxLength)+' 'str(boxLength)
        else:
            assert('Please ensure the simulation cell is set to be either tailored, or cubic.')
        
        changeBoxSettings(ufl,diag)
        
        MG(toolsDir,ufl)
        
        DD2OvitoVTK(toolsDir,ufl)
        
        copyConfigAndResults(ufl,dataLoc)
        
    #run the CRSS detection
    
    #######################CRSS DETECTION & CONVERGENCE STUDY##################
    CRSS_known=1; #switch to 1 to disable CRSS convergence
    if CRSS_known==1:
        CRSS_set=150e6# [Pa]
    startStress_Pa=120e6 #First stress to attempt for CRSS
    endStress_Pa=500e6 #last stress to attempt for CRSS
    origin=[0,0,0]
    ssID=0
    numNodesMat=[1,5,10,20,30,40,50,60]
    #    numNodesMat=[1,5]
    exitFaceID=2
    initialDipoleGlideSteps=1
    dipoleHeight=14
    partialsActive=0
    numSimSteps=1000
    PBCs='8 8 8'
    ewaldLengthFactor='8'
    temperature=300
    langevinThermostat=0
    langevinSeed=0
    quadPerLengthMat=[0.1,0.3,0.6,0.9,1]
    #    quadPerLengthMat=[1]
    numRelaxSteps=10


    if CRSS_known==0:
        runCRSS_convergence(ufl,startStress_Pa,endStress_Pa,origin,ssID,numNodesMat,exitFaceID,initialDipoleGlideSteps,dipoleHeight,partialsActive,numSimSteps,PBCs,ewaldLengthFactor,temperature,langevinThermostat,langevinSeed,quadPerLengthMat,numRelaxSteps,dataLoc)
    
    
    
    #store the data




    return()


if __name__ == "__main__":
    main()
