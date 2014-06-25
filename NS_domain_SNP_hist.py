#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np

bed_file = sys.argv[1]


class bedDetail:
    def __init__(self, GStart, GStop, EStart, EStop, GeneName, ExonName):
        self.GS = GStart
        self.GP = GStop
        self.ES = EStart
        self.EP = EStop
        self.GN = GeneName
        self.EN = ExonName



beds = []
for line in open(bed_file, 'r'):
    l = line.split("\t")
    gstart = int(l[1])
    gstop = int(l[2])
    estart = int(l[-2])
    estop = int(l[-1])
    feature = str(l[3])
    f = feature.split(" ")
    gene = str(f[0])
    d = f[1:]
    domain = str(' '.join(d))
    x = bedDetail(gstart, gstop, estart, estop, gene, domain)
    beds.append(x)
    

nDomains = len(beds)
#print nDomains
results_file = sys.argv[2]

names = []
results = []
for line in open(results_file, "rU"):
    line_results = line.split("\t")
    resultsL = []
    if line.startswith('H'):
        for l in line_results:
            names.append(l.replace("\n",""))        
    else:
        for l in line_results:
            resultsL.append(int(l))
        results.append(resultsL)
    
#print results
#print names
animals = len(results[0]) # calculate number of animals
ind = np.arange(nDomains)    # the x locations for the groups
width = 0.90       # the width of the bars: can also be len(x) sequence
ax = plt.subplot2grid((24,12), (0,0), rowspan=18, colspan=11)
score50 = []

for f in results:
    score50.append(0)

    
#colors = ["#FF0000", "#009999", "#9FEE00", "#BF3030", "#1D7373", "#86B32D", "#A60000",\
#"#006363", "#679B00", "#FF4040", "#33CCCC", "#B9F73E", "#FF7373", "#5CCCCC", "#C9F76F"]


colors = ["k", \
"#a6cee3", \
"#1f78b4", \
"0.8", \
"#b2df8a", \
"#33a02c", \
"#fb9a99", \
"0.6", \
"#e31a1c", \
"#fdbf6f", \
"#ff7f00", \
"#cab2d6", \
"#6a3d9a", \
"#ffff99", \
"#a6cee3", \
"#1f78b4", \
"0.8", \
"#b2df8a", \
"#33a02c", \
"#fb9a99", \
"0.6", \
"#e31a1c", \
"#fdbf6f", \
"#ff7f00", \
"#cab2d6", \
"#6a3d9a", \
"#ffff99", \
"#b15928"]
				

for i in range(animals):
    pi = "p" + str(i)
    score75 = []
    for f in results:
        score75.append(f[i])
    #print pi, score75
    pi = ax.bar(ind, score75, width, bottom=score50,color=colors[i], label=names[i])
    lists = []
    lists.append(score75)
    lists.append(score50)
    score50 = map(sum, zip(*lists))

plt.ylabel('Number of SNPs per domain')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':11})
ax.tick_params(axis='x', which='major', labelsize=8)
ax.tick_params(axis='x', which='minor', labelsize=6)
plt.xlim([0,len(score75)])
domains = []
for b in beds:
    domains.append(str(b.EN))
plt.xticks(ind+width/2., (domains), rotation=90)
plt.show()
#plt.savefig(savename, format='pdf')

