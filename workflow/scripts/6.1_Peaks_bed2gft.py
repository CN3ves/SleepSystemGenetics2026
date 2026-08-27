# Script to split all peaks into smaller regions and convers bed to gft

import sys
import math

FILE = sys.argv[1]
OUTPUT = sys.argv[2]

original_size =[]
new_size=[]

Data = {}
PeakID = 0
# Collect all feature data into a dictionary
with open(FILE, "r") as input:
  for line in input:
    line = line.replace("\n","").split("\t")
    Chr = line[0]
    Start = int(line[1])+1 # shift base gtf are 1-indexed files while bed are 0-indexed
    End = int(line[2])+1
    
    PeakID += 1
    
    original_size.append(End - Start)
    # Split peaks longer than 250 into equally sized regions
    if original_size[-1] > 250: 
      split_peak = [i for i in  range(Start, End, 250)]
      
      if (split_peak[-1] - split_peak[-2]) < 100: # if last region is very small, split the last two equally
        step = math.floor((split_peak[-1] - split_peak[-2]) / 2)
        split_peak[-1] = split_peak[-2] + step
      
      split_peak.append(End)
      
    else:
      split_peak = [Start, End]
        
    for i in range(1,len(split_peak)):
        new_size.append(split_peak[i] - split_peak[i-1])
        if Chr not in Data: 
          Data[Chr] = {}
        #format: [seqname, source, feature, start, end, score, strand, frame, attribute]
        ID=f"{PeakID}.{i}"
        Data[Chr][ID]=[Chr, "Genrich", "Peak", str(split_peak[i-1]), str(split_peak[i]), ".", "+", ".", f"Peak_ID {Chr}_{ID};"]

print(f"Last Peak: {PeakID}")    

print(f"Original range: {min(original_size)} - {max(original_size)}\nSplit range: {min(new_size)} - {max(new_size)}")
# Print modified feature data
with open(OUTPUT, "w") as output:
  chrs = sorted(Data.keys())
  for Chr in chrs:
    peaks = sorted(Data[Chr].keys())
    for Peak in peaks:
      output.write("\t".join(Data[Chr][Peak])+"\n")