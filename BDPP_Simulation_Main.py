
import SimFunctions
import SimRNG 
import SimClasses
import pandas as pd
import numpy as np
import statistics as stats
import math

ZSimRNG = SimRNG.InitializeRNSeed()

# %% Simulation settings ####################################################

InocQueue = SimClasses.FIFOQueue()  # queue Incolation
MainQueue = SimClasses.FIFOQueue()  # queue forMain Vessel
CentQueue = SimClasses.FIFOQueue()  # queue for centrifuge
ChromQueue = SimClasses.FIFOQueue()  # queue for chromatography
FilQueue = SimClasses.FIFOQueue()  # queue for filtration

Wait = SimClasses.DTStat()

InocServer = SimClasses.Resource()  # server for Inoculum Tanks
MainServer = SimClasses.Resource()  # server for Main Vessel
CentServer = SimClasses.Resource()  # server for Centrifuge Ewuipment
ChromServer = SimClasses.Resource()  # server for Chromatography Equipment
FilServer = SimClasses.Resource()  # server for Filtration Equipment
Calendar = SimClasses.EventCalendar()

TheCTStats = []
TheDTStats = []
TheQueues = []
TheResources = []

TheDTStats.append(Wait)

TheQueues.append(InocQueue)
TheQueues.append(MainQueue)
TheQueues.append(CentQueue)
TheQueues.append(ChromQueue)
TheQueues.append(FilQueue)

TheResources.append(InocServer)
TheResources.append(MainServer)
TheResources.append(CentServer)
TheResources.append(ChromServer)
TheResources.append(FilServer)

InocServer.SetUnits(8)
MainServer.SetUnits(8)
CentServer.SetUnits(4)
ChromServer.SetUnits(4)
FilServer.SetUnits(4)

alpha = 1.5 #used to calculate impurity level at main fermentation level
Beta = 1.0 #used to calculate variance of impurity at fermentation level
Beta_ = [0.5, 1.0, 1.5]
gamma = 0.15
# WarmUp = 100
# MeanTBA = 1
# MeanST = 5.0*1.05
InocST = 24
MainST = 72
CentST = 2.5
ChromST = 2.0
FilST = 2.0

rep = 100
breps = 200
MaxBatches = 200
RunLength = 5000

# AllWaitMean = []
# MainQueueMean = []
# CentQueueMean = []
# ChromQueueMean = []
# FilQueueMean = []
# AllMainQueueNum = []
# AllCentQueueNum = []
# AllChromQueueNum = []
# AllFilQueueNum = []
# AllMainServerMean = []
# AllCentServerMean = []
# AllChromServerMean = []
# AllFilServerMean = []


#%% Bootstraping #########################################################

df_ferm = pd.read_excel("C:/Users/prami/Downloads/Data BioManufacturing _ Project 1 _ Pramit Gopal Yeole _ Andrew Le _ IE 7215.xlsx", sheet_name="Fermentation")
df_chrom = pd.read_excel("C:/Users/prami/Downloads/Data BioManufacturing _ Project 1 _ Pramit Gopal Yeole _ Andrew Le _ IE 7215.xlsx", sheet_name="Chromatography")

# =============================================================================
# Distribution obtained from @Risk for
# Initial Biomass -- X0 ~ N(4.8878935, 1.599819**2)
#
# Main Fermentation
# (uT + Ep) ~ N(1.7931284, 0.2686605**2)
# Ei ~ N(0, 0.2686605**2)
# Xf = X0*exp(uT + Ep)
# If = Xf*alpha*exp(Ei)

# Centrifuge
# Q ~ Unif(0.4, 0.5) #constant
# Xc = Xf

# Chromatography
# Xp = Qp*Xc | Qp ~ Unif(0.6571625, 0.9178249)
# Ip = Qi*Ic | Qi ~ Unif(0.2424452, 0.6645715)

# Filtration
# Xfr = Xp
# Ifr = Qfr*Ip | Qfr ~ Unif(0.99, 1)
# =============================================================================

# data for sampling with replacement
X0_data = df_ferm['Initial Biomass X_0'].values
Xf_data = df_ferm['Protein X_F'].values
If_data = df_ferm['Impurity I_F'].values

MeanProteinlevel = []
MeanImpurity = []
MeanCycleTime = []
StdProteinlevel = []
StdImpurity = []
StdCycleTime = []
Beta_Mean = {}
Beta_Std = {}
Beta_MeanCI = {}
Beta_StdCI = {}
Beta_error_Mean = {}
Beta_error_Std = {}

#%% Simulation Process ####################################################

def Arrival(i):
    SimFunctions.Schedule(Calendar, "Arrival", SimRNG.Expon(2, 24)) #unif dist is the random arrival process
    Sample = SimClasses.Entity() # only one entity created for the batch i.e. Antigen A
    Sample.SampleNo = i + 1
    Sample.SampleMass = SimRNG.Normal(X0_mean, X0_var, 1)
    Biomass_X0.append((Sample.SampleNo, Sample.SampleMass)) #appending data of the batch and the biomass generated
    InocQueue.Add(Sample) #adding the sample data to be the Inoculum Ferm queue for subsequent processing
    if InocServer.Busy < InocServer.NumberOfUnits: # check whether can get immediate service
        InocQueue.Remove()  # when start service remove from queue
        InocServer.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndOfInocFermentation", InocST, Sample)  # schedule the end of inoculum fermentation

def EndOfInocFermentation(DepartingSample):
    # Wait.Record(SimClasses.Clock - DepartingSample.CreateTime)
    if InocQueue.NumQueue() > 0:
        Sample = InocQueue.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndOfInocFermentation", InocST, Sample)
        SimFunctions.SchedulePlus(Calendar, "ArrivalMainFermentation", 0, DepartingSample)
    else:
        InocServer.Free(1)
        SimFunctions.SchedulePlus(Calendar, "ArrivalMainFermentation", 0, DepartingSample)

def ArrivalMainFermentation(ArrivingSample):
    MainQueue.Add(ArrivingSample)
    if MainServer.Busy < MainServer.NumberOfUnits: # check whether can get immediate service
        MainQueue.Remove()  # when start service remove from queue
        MainServer.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndOfMainFermentation", MainST, ArrivingSample)

def EndOfMainFermentation(DepartingSample):
    # Wait.Record(SimClasses.Clock - DepartingSample.CreateTime)
    if MainQueue.NumQueue() > 0:
        Sample = MainQueue.Remove()
        X_f = DepartingSample.SampleMass*(np.exp(SimRNG.Normal(exp_Xf_mean, exp_Xf_var, 3)))
        MainFermXF.append(X_f)
        I_f = X_f*alpha*(np.exp(SimRNG.Normal(exp_If_mean, exp_If_var, 3)))
        MainFermIF.append(I_f)
        SimFunctions.SchedulePlus(Calendar, "EndOfMainFermentation", MainST, Sample)
        SimFunctions.SchedulePlus(Calendar, "ArrivalCentrifuge", 0, DepartingSample)
    else:
        MainServer.Free(1)
        SimFunctions.SchedulePlus(Calendar, "ArrivalCentrifuge", 0, DepartingSample)

def ArrivalCentrifuge(ArrivingSample):
    CentQueue.Add(ArrivingSample)
    if CentServer.Busy < CentServer.NumberOfUnits: # check whether can get immediate service
        CentQueue.Remove()  # when start service remove from queue
        CentServer.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndOfCentrifuge", CentST, ArrivingSample)

def EndOfCentrifuge(DepartingSample):
    # Wait.Record(SimClasses.Clock - DepartingSample.CreateTime)
    I_c = SimRNG.Uniform(0.4, 0.5, 4)*MainFermIF[DepartingSample.SampleNo-1]
    CentIC.append(I_c)
    if CentQueue.NumQueue() > 0:
        Sample = CentQueue.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndOfCentrifuge", CentST, Sample)
        SimFunctions.SchedulePlus(Calendar, "ArrivalChromatography", 0, DepartingSample)
    else:
        CentServer.Free(1)
        SimFunctions.SchedulePlus(Calendar, "ArrivalChromatography", 0, DepartingSample)

def ArrivalChromatography(ArrivingSample):
    ChromQueue.Add(ArrivingSample)
    if ChromServer.Busy < ChromServer.NumberOfUnits: # check whether can get immediate service
        ChromQueue.Remove()  # when start service remove from queue
        ChromServer.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndOfChromatography", ChromST, ArrivingSample)

def EndOfChromatography(DepartingSample):
    # Wait.Record(SimClasses.Clock - DepartingSample.CreateTime)
    X_p = SimRNG.Uniform(0.6571625, 0.9178249, 6)*MainFermIF[DepartingSample.SampleNo-1]
    ChromXP.append(X_p)
    I_p = SimRNG.Uniform(0.2424452, 0.6645715, 7)*CentIC[DepartingSample.SampleNo-1]
    ChromIP.append(I_p)
    if ChromQueue.NumQueue() > 0:
        Sample = ChromQueue.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndOfChromatography", ChromST, Sample)
        SimFunctions.SchedulePlus(Calendar, "ArrivalFiltration", 0, DepartingSample)
    else:
        ChromServer.Free(1)
        SimFunctions.SchedulePlus(Calendar, "ArrivalFiltration", 0, DepartingSample)

def ArrivalFiltration(ArrivingSample):
    FilQueue.Add(ArrivingSample)
    if FilServer.Busy < FilServer.NumberOfUnits: # check whether can get immediate service
        FilQueue.Remove()  # when start service remove from queue
        FilServer.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndOfFiltration", ChromST, ArrivingSample)
        
def EndOfFiltration(DepartingSample):
    # Wait.Record(SimClasses.Clock - DepartingSample.CreateTime)
    I_Fr = SimRNG.Uniform(0.99, 1, 10) * ChromIP[DepartingSample.SampleNo-1]
    FilIFr.append(I_Fr)
    CycleTime.append(SimClasses.Clock - DepartingSample.CreateTime)
    if FilQueue.NumQueue() > 0:
        Sample = FilQueue.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndOfFiltration", FilST, Sample)
    else:
        FilServer.Free(1)
        
#%% Main ##################################################################

# sensitivity analysis with varying Beta values
for b in Beta_:
    for reps in range(breps): #bootstraping 200 samples
        X0_sample = np.random.choice(X0_data, size=len(X0_data)) #array
        Xf_sample = np.random.choice(Xf_data, size=len(Xf_data))
        
        Biomass_X0 = []
        MainFermXF = []
        MainFermIF = []
        CentIC = []
        ChromXP = []
        ChromIP = []
        FilIFr = []
        CycleTime = []
        Inoc_CycleTime = []
        Main_CycleTime = []
        Cent_CycleTime = []
        Chrom_CycleTime = []
        Fil_CycleTime = []
        
        #running the no. of reps = 100
        for j in range(rep):
            #get the distribution parameters
            X0_mean = np.mean(X0_sample)
            X0_var = np.var(X0_sample)
            
            #to get parameters for the exp for Xf
            exp_Xf_sample = np.log(Xf_sample) - np.log(X0_sample)
            exp_Xf_mean = np.mean(exp_Xf_sample)
            exp_Xf_var = np.var(exp_Xf_sample)
            exp_If_mean = 0
            exp_If_var = (b**2) * exp_Xf_var
            
            SimFunctions.SimFunctionsInit(Calendar, TheQueues, TheCTStats, TheDTStats, TheResources)
            SimFunctions.Schedule(Calendar, "Arrival", SimRNG.Expon(2, 24)) #this is the arrival
            SimFunctions.Schedule(Calendar, "EndSimulation", RunLength)
            i = 0
            
            NextEvent = Calendar.Remove()
            SimClasses.Clock = NextEvent.EventTime
            if NextEvent.EventType == "Arrival":
                Arrival(i)
                i=i+1
            elif NextEvent.EventType == "EndOfInocFermentation":
                depart_one = NextEvent.WhichObject
                EndOfInocFermentation(depart_one)
            elif NextEvent.EventType == "ArrivalMainFermentation":
                depart_one = NextEvent.WhichObject
                ArrivalMainFermentation(depart_one)
            elif NextEvent.EventType == "EndOfMainFermentation":
                depart_one = NextEvent.WhichObject
                EndOfMainFermentation(depart_one)
            elif NextEvent.EventType == "ArrivalCentrifuge":
                depart_one = NextEvent.WhichObject
                ArrivalCentrifuge(depart_one)
            elif NextEvent.EventType == "EndOfCentrifuge":
                depart_one = NextEvent.WhichObject
                EndOfCentrifuge(depart_one)
            elif NextEvent.EventType == "ArrivalChromatography":
                depart_one = NextEvent.WhichObject
                ArrivalChromatography(depart_one)
            elif NextEvent.EventType == "EndOfChromatography":
                depart_one = NextEvent.WhichObject
                EndOfChromatography(depart_one)
            elif NextEvent.EventType == "ArrivalFiltration":
                depart_one = NextEvent.WhichObject
                ArrivalFiltration(depart_one)
            elif NextEvent.EventType == "EndOfFiltration":
                depart_one = NextEvent.WhichObject
                EndOfFiltration(depart_one)
            
            while NextEvent.EventType != "EndSimulation":
                if len(FilIFr) >= MaxBatches:
                    break
                NextEvent = Calendar.Remove()
                SimClasses.Clock = NextEvent.EventTime
                if NextEvent.EventType == "Arrival":
                    Arrival(i)  
                    i=i+1
                elif NextEvent.EventType == "EndOfInocFermentation":
                    depart_one = NextEvent.WhichObject
                    EndOfInocFermentation(depart_one)
                elif NextEvent.EventType == "ArrivalMainFermentation":
                    depart_one = NextEvent.WhichObject
                    ArrivalMainFermentation(depart_one)
                elif NextEvent.EventType == "EndOfMainFermentation":
                    depart_one = NextEvent.WhichObject
                    EndOfMainFermentation(depart_one)
                elif NextEvent.EventType == "ArrivalCentrifuge":
                    depart_one = NextEvent.WhichObject
                    ArrivalCentrifuge(depart_one)
                elif NextEvent.EventType == "EndOfCentrifuge":
                    depart_one = NextEvent.WhichObject
                    EndOfCentrifuge(depart_one)
                elif NextEvent.EventType == "ArrivalChromatography":
                    depart_one = NextEvent.WhichObject
                    ArrivalChromatography(depart_one)
                elif NextEvent.EventType == "EndOfChromatography":
                    depart_one = NextEvent.WhichObject
                    EndOfChromatography(depart_one)
                elif NextEvent.EventType == "ArrivalFiltration":
                    depart_one = NextEvent.WhichObject
                    ArrivalFiltration(depart_one)
                elif NextEvent.EventType == "EndOfFiltration":
                    depart_one = NextEvent.WhichObject
                    EndOfFiltration(depart_one)
            # print(CycleTime)    
        MeanProteinlevel.append(np.mean(ChromXP)) #calculating mean
        MeanImpurity.append(np.mean(np.asarray(FilIFr)))
        MeanCycleTime.append(np.mean(np.asarray(CycleTime)))
        
        StdProteinlevel.append(np.std(ChromXP))
        StdImpurity.append(np.std(FilIFr))
        StdCycleTime.append(np.std(CycleTime))
    
        
    
    # used to calculate the values for a single Beta value
    Boot_MeanProtein = np.mean(MeanProteinlevel)
    Boot_MeanImpurity = np.mean(MeanImpurity)
    Boot_MeanCycleTime = np.mean(MeanCycleTime)
    
    Boot_stdProtein = np.mean(StdProteinlevel)
    Boot_stdImpurity = np.mean(StdImpurity)
    Boot_stdCycleTime = np.mean(StdCycleTime)
    
    # respective errors
    boot_error_meanprotein = 1.96*np.std(MeanProteinlevel)/math.sqrt(len(MeanProteinlevel))
    boot_error_meanimpurity = (1.96*np.std(MeanImpurity)/math.sqrt(len(MeanImpurity)))
    boot_error_meancycletime = 1.96*np.std(MeanCycleTime)/math.sqrt(len(MeanCycleTime))
    
    boot_error_stdprotein = 1.96*np.std(StdProteinlevel)/math.sqrt(len(StdProteinlevel))
    boot_error_stdimpurity = (1.96*np.std(StdImpurity)/math.sqrt(len(StdImpurity)))
    boot_error_stdcycletime = 1.96*np.std(StdCycleTime)/math.sqrt(len(StdCycleTime))

    # Calculating the 95% confidence interval for 200 bootstraped samples
    MeanProteinCI = [Boot_MeanProtein - 1.96*np.std(MeanProteinlevel)/math.sqrt(len(MeanProteinlevel)), Boot_MeanProtein + 1.96*np.std(MeanProteinlevel)/math.sqrt(len(MeanProteinlevel))]
    MeanImpurityCI = [Boot_MeanImpurity - (1.96*np.std(MeanImpurity)/math.sqrt(len(MeanImpurity))), Boot_MeanImpurity + 1.96*np.std(MeanImpurity)/math.sqrt(len(MeanImpurity))]
    MeanCycleTimeCI = [Boot_MeanCycleTime - 1.96*np.std(MeanCycleTime)/math.sqrt(len(MeanCycleTime)), Boot_MeanCycleTime + 1.96*np.std(MeanCycleTime)/math.sqrt(len(MeanCycleTime))]
    
    StdProteinCI = [Boot_stdProtein - 1.96*np.std(StdProteinlevel)/math.sqrt(len(StdProteinlevel)), Boot_stdProtein + 1.96*np.std(StdProteinlevel)/math.sqrt(len(StdProteinlevel))]
    StdImpurityCI = [Boot_stdImpurity - (1.96*np.std(StdImpurity)/math.sqrt(len(StdImpurity))), Boot_stdImpurity + 1.96*np.std(StdImpurity)/math.sqrt(len(StdImpurity))]
    StdCycleTimeCI = [Boot_stdCycleTime - 1.96*np.std(StdCycleTime)/math.sqrt(len(StdCycleTime)), Boot_stdCycleTime + 1.96*np.std(StdCycleTime)/math.sqrt(len(StdCycleTime))]
    
    #Calculating the Variance of the paarmeters
    VarMeanProteinLevel = np.var(MeanProteinlevel)
    VarMeanImpurity = np.var(MeanImpurity)
    VarMeanCycleTime = np.var(MeanCycleTime)
    
    VarStdProteinLevel = np.var(StdProteinlevel)
    VarStdImpurity = np.var(StdImpurity)
    VarStdCycleTime = np.var(StdCycleTime)
    
    #appended as Protein, Impurity, CycleTime
    Beta_Mean[b]= [Boot_MeanProtein, Boot_MeanImpurity, Boot_MeanCycleTime]
    Beta_Std[b] = [Boot_stdProtein, Boot_stdImpurity, Boot_stdCycleTime]
    
    Beta_MeanCI[b] = [MeanProteinCI, MeanImpurityCI, MeanCycleTimeCI]
    Beta_StdCI[b] = [MeanProteinCI, MeanImpurityCI, MeanCycleTimeCI]
    
    Beta_error_Mean[b] = [boot_error_meanprotein, boot_error_meanimpurity, boot_error_meancycletime]
    Beta_error_Std[b] = [boot_error_stdprotein, boot_error_stdimpurity, boot_error_stdcycletime]
  