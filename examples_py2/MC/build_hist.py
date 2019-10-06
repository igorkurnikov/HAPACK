# The script reads the tortion angle trajectory file ene.out and builds the histogram
# of the angle values to calculate probability distribution
# output is to the file hist.out
# probability distribution is not normalized
fmin = 0.0 #minimum angle value
fmax = 360.0 #maximum angle value
nbin = 20 # number of bins in the histogram
delt = (fmax - fmin)/nbin #bin with
prob = [] #an array to collect angle values to calculate probability
npr_freq  = 100 # frequency to print info to console

for i in range(nbin):
  prob.append(0)
finp = open("ene.out","r") #read MC trajectory file
str = finp.readline()
istep = 0
while( len(str)  > 0 ): #cycle through the whole traj file
  words = str.split()
  ang = float(words[1])
  istep= istep+1
  if( istep % npr_freq == 0):
    print istep, "  ",ang
  idx =  int( (ang - fmin)/delt) #calculate  to which bin ang belongs 
  prob[idx] = prob[idx] + 1 #add a count to that bin
  str = finp.readline()
finp.close()
fout = open("hist.out","w")
for i in range(nbin):
  ff = fmin + delt*( i + 0.5) #add 0.5 to plot hist value in the middle of a bin
  print >>fout, i, ff, prob[i]
fout.close()
  
  