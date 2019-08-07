import random
import math
import sys
from bisect import bisect_left

"""
Writen by Michael Bottery 2015
"""

class gbStack:
    def __init__(self, totalTime, seed):
        self.heaptimes = []
        self.heapreactions = []
        self.heapProps = []
        self.totalTime = totalTime
        self.lastReaction = None
        self.heaptimes.append(float("inf"))
        self.output = False
        random.seed(seed)

    
    def checkElement(self, reaction):
        return self.heapreactions.index(reaction)

        
    def countElemets(self):
        return len(self.heapreactions)

    
    def getProp(self, reaction):
        position = self.heapreactions.index(reaction)
        return self.heapProps[position]

        
    def returnElement(self, i):
        return self.heapreactions[i]

    
    def addElement(self, time, reaction, prop = -1.0):
        #if self.output: print("Adding {0} of {1} to the list at time {2} with prop {3}".format(
        #                       time, reaction, self.totalTime, prop)),
		
        Time = time + self.totalTime
		
        pos = bisect_left(self.heaptimes, Time)
        
        self.heaptimes.insert(pos,Time);
        self.heapreactions.insert(pos,reaction);
        self.heapProps.insert(pos,prop);

        #if self.output: print(" -- it went in at posn {0}".format(pos))


    def updateElement(self, reaction, oldProp, newProp):
        position = self.heapreactions.index(reaction)
        tau = self.heaptimes[position]
        
        if math.fabs(oldProp - self.heapProps[position]) > 0.0000001:
            print("OldProp: {0}, stack prop = {1}, id: {2}".format(oldProp,
                   self.heapProps[position], reaction))
            sys.exit("**************** Sync lost ****************")
        
        if(oldProp == newProp):
            pass
            #if self.output: 
            #    print("OldProp: {0}, stack prop = {1}, id: {2}, DO NOTHING".format(oldProp,
            #           self.heapProps[position], reaction))
        else:
            del self.heaptimes[position]
            del self.heapreactions[position]
            del self.heapProps[position]
            
            #if self.output:
            #    print("OldProp {0} {1} {2} {3}").format(oldProp, newProp, 
            #           (tau-self.totalTime), (oldProp/newProp))
            #    print("{0}: Updating {1}: old time was {2}".format(self.lastReaction,
            #           reaction, tau)),
            if tau == float("inf"):
                tau = -1.0*math.log(random.random())/newProp
                #if self.output:
                #    print(", we generated a new one and"),
            else:
                while tau == self.totalTime:
                    tau += 0.1
                    sys.exit("**************** Time clash ****************")
                    #if self.output:
                     #   print(", CORRECTED to {0}".format(tau)),
                #if self.output:
                 #   print(", we rescaled the old one and"),
                tau = (tau-self.totalTime)*(oldProp/newProp)
            #if self.output:
              #  print(", new time is {0}".format(tau + self.totalTime))
            self.addElement(tau, reaction, newProp)

    
    def removeOne(self, type):
        count = 0
        returnTime = -1.0
        for i in xrange(len(self.heapreactions)):
            if self.heapreactions[i] == type:
                count+=1
        
        if count > 0:
            q = 1 + math.floor(random.random()*count) #floor(x), largest int <= x
            count = 0 
            i = 0
            while count < q:
                if self.heapreactions[i] == type:
                    count+=1
                i+=1
            i-=1

            returnTime = self.heaptimes[i]
            
            del self.heapreactions[i]
            del self.heaptimes[i]
            del self.heapProps[i]
            
        return returnTime


    def dump(self):
        print("position\theaptime\theapreaction\theapProp")
        for i in xrange(len(self.heapreactions)):
            print("{0}\t\t{1:.4f}\t\t{2}\t\t{3:.4f}".format(i, self.heaptimes[i],
                   self.heapreactions[i], self.heapProps[i]))


    def nextReaction(self):
        self.lastReaction = self.heapreactions[0]
        if self.output:
            print("Chosen reaction: {0}".format(self.lastReaction))
        del self.heapreactions[0]
        self.totalTime = self.heaptimes[0]
        del self.heaptimes[0]
        del self.heapProps[0]
        
        return self.lastReaction


    def readNextTime(self):
        return self.heaptimes[0]


    def tnextGaussian(self):
        x1 = 1.0
        x1 = 1.0
        w = x1*x1+x2*x2
        
        while x > 1:
            x1 = 2*random.random()-1.0
            x2 = 2*random.random()-1.0
            w = x1*x1+x2*x2
        
        return math.sqrt(-2.0*mathlog(w)/w)*x1


    def timeShift(self, shift):
        self.totalTime -= shift
        
        for i in xrange(len(self.heapreactions)):
            r = self.heapreactions[i]
            t = self.heaptimes[i]
            p = self.heapProps[i]
            del self.heapreactions[i]
            del self.heaptimes[i]


    def clean(self):
        self.heapreactions = []
        self.heaptimes = []
        self.heapProps = []
        self.heaptimes.append(float("inf"))
