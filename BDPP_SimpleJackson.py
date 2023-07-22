import SimFunctions
import SimRNG
import SimClasses
import pandas as pd
import numpy as np
import statistics as stats
import math
import scipy.stats as sp_stats


# Initialization
ZSimRNG = SimRNG.InitializeRNSeed()

sample_times = {}

# Create Nodes for Jackson Network
class Node:
    def __init__(self, name, service_time, stream, num_units=1):
        self.name = name
        self.service_time = SimRNG.Expon(service_time, stream)
        self.queue = SimClasses.FIFOQueue()
        self.resource = SimClasses.Resource()
        self.resource.SetUnits(num_units)

nodes = [
    Node("InocFermentation", 24, 2, 4),
    Node("MainFermentation", 72, 4, 8),
    Node("Centrifuge", 2.5, 6, 8),
    Node("Chromatography", 2.0, 8, 8),
    Node("Filtration", 2.0, 10, 8)
]

# Function to execute a stage in the process
def process_stage(stage, sample):
    # Perform the logic for each stage
    if stage == 0:
        # InocFermentation
        pass
    elif stage == 1:
        # MainFermentation
        pass
    elif stage == 2:
        # Centrifuge
        pass
    elif stage == 3:
        # Chromatography
        pass
    elif stage == 4:
        # Filtration
        pass

    # Check if there's a next stage
    if stage + 1 < len(nodes):
        next_node = nodes[stage + 1]
        next_node.queue.Add(sample)
        if next_node.resource.Busy < next_node.resource.NumberOfUnits:
            next_node.queue.Remove()
            next_node.resource.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndOfStage", next_node.service_time, (stage + 1, sample))
    else:
        sample_times[sample]['end_time'] = SimClasses.Clock

# Schedule events
def arrival(stage, sample):
    if stage == 0:
        sample_times[sample] = {'start_time': SimClasses.Clock}
    SimFunctions.SchedulePlus(Calendar, "EndOfStage", SimRNG.Expon(1,24), (stage, sample))

def end_of_stage(stage, sample):
    process_stage(stage, sample)

    # Process the rest of the stages
    current_node = nodes[stage]
    if current_node.queue.NumQueue() > 0:
        next_sample = current_node.queue.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndOfStage", current_node.service_time, (stage, next_sample))
    else:
        current_node.resource.Free(1)

# Main simulation loop
for _ in range(200):
    # Initialize Calendar and other required instances
    Calendar = SimClasses.EventCalendar()
    TheQueues = []
    TheCTStats = []
    TheDTStats = []
    TheResources = []
    for node in nodes:
        TheQueues.append(node.queue)
        TheCTStats.append(SimClasses.CTStat())
        TheDTStats.append(SimClasses.DTStat())
        TheResources.append(node.resource)

    # Simulate the process for each sample
    for i in range(200):
        arrival(0, i)

    # Run the simulation
    while Calendar.N() > 0:
        NextEvent = Calendar.Remove()
        SimClasses.Clock = NextEvent.EventTime
        stage, sample = NextEvent.WhichObject
        end_of_stage(stage, sample)
        
cycle_times = []
for sample in sample_times:
    start_time = sample_times[sample]['start_time']
    end_time = sample_times[sample]['end_time']
    cycle_time = end_time - start_time
    cycle_times.append(cycle_time)

average_cycle_time = sum(cycle_times) / len(cycle_times)
print("Average Cycle Time:", average_cycle_time)

# Calculate mean and standard deviation
mean_cycle_time = np.mean(cycle_times)
std_cycle_time = np.std(cycle_times, ddof=1)

# Calculate the confidence interval
confidence_level = 0.95
degrees_of_freedom = len(cycle_times) - 1
t_value = sp_stats.t.ppf((1 + confidence_level) / 2, degrees_of_freedom)

margin_of_error = t_value * (std_cycle_time / math.sqrt(len(cycle_times)))
confidence_interval = (mean_cycle_time - margin_of_error, mean_cycle_time + margin_of_error)

print("Confidence Interval (95%):", confidence_interval)

