from functions.psuff import * 

data = readSegyFile("data/seismic.segy") 

showBinaryHeader(data)
showTraceHeader(data)

perc = 99
ishot = 10
fshot = 40

shots = buildShotWindow(data,ishot,fshot)
plotShotWindow(data,shots,perc,"figures/shots.png")

plotAcquisitionGeometry(data,ishot,fshot,"figures/geometria.png")

poffset = showOffsetPossibilities(data)
section = buildOffsetSection(data,poffset[0])
plotOffsetSection(data,section,poffset[0],perc,"figures/offset.png")

line = showCMPsPossibilities(data)
plotCMPsPosition(data,line,"figures/cmps.png")

gatherCMP = buildCMPSection(data,line[100])
plotCMPGather(data,gatherCMP,line[100],perc,"figures/gatherCMP.png")

