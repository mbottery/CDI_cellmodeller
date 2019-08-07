import numpy as np
import sys
import os
sys.path.append('.')
import cPickle
import glob
import operator

class radial():
    # Class that collects radial data from pickle files
    def __init__(self, data):
        self.states = data['cellStates']
        self.radialPosition = {0:[],1:[],2:[],3:[]} 
        for (id,s) in self.states.iteritems():
            self.radialPosition[s.cellType].append(np.sqrt(s.pos[0]**2+s.pos[1]**2))
            
    def calc_radius(self):
        return max(self.radialPosition[0]+self.radialPosition[1]+self.radialPosition[2]+self.radialPosition[3])
    
    def radial_fraction(self,max_r = 400,bin_size = 20):        
        r = range(0,max_r,bin_size)
        hist_celltype_0 = np.histogram(self.radialPosition[0],bins = r)[0].astype(float)
        hist_celltype_1 = np.histogram(self.radialPosition[1],bins = r)[0].astype(float)
        hist_celltype_2 = np.histogram(self.radialPosition[2],bins = r)[0].astype(float)
        hist_celltype_3 = np.histogram(self.radialPosition[3],bins = r)[0].astype(float)
        ratio_cellType_0 = hist_celltype_0/(hist_celltype_0+hist_celltype_1+hist_celltype_2+hist_celltype_3)
        return r[:-1],ratio_cellType_0
        
        
def openPickle(pickle):
    print('opening: '+ pickle)
    data = False
    try:
        data = cPickle.load(open(pickle,'r'))
    except:
        print 'Failed to open pickle file'
        print sys.exc_info()[0]
    return data

def main():
    # open the final step file (step-00700.pickle) in each of the directories
    # calculate radius of the colonies
    # bin the fraction of cells based on radial position
    # return the data to a main data table with all other replicates
    # table must have all the parameters, machine, replicate number, fraction
    # cell type 1, fraction cell type 2 and radius of outer edge of 
    # annulus.
    # write out the table so I can plot it in R
    indirs = sys.argv[1:]
    
    radius_table = open('growth_difference_radius.txt','w')
    radius_table.write('growth_difference,density,rep,radius\n')
    
    radial_fraction_table = open('growth_difference_radial_fraction.txt','w')
    radial_fraction_table.write('growth_inhibited,density,rep,radius,ratio_0\n')
    
    growth_difference = {'0':'0','1':'0.01','2':'0.05','3':'0.1','4':'0.2'}
    density = {'0':'20000','1':'2000','2':'200'}
    
    codes=[]
    
    for d in range(len(indirs)): # machine
        for dirs in os.listdir(indirs[d]):  # simulations 
            print os.path.join(indirs[d],dirs)
            code = dirs[18:21]
            print 'simulation: {0}'.format(code)
            codes.append(code)
            gd,de = code[0],code[2]
            rep = indirs[d][3:]
            print 'rep: {0}'.format(rep)
            for filename in os.listdir(os.path.join(indirs[d],dirs)): # pickle steps
                if filename == "step-00700.pickle":
                    data = openPickle(os.path.join(indirs[d],dirs,filename))
                    if data:
                        r = radial(data)
                        radius = r.calc_radius()
                        bins,fraction_cellType_0 = r.radial_fraction()
                    else:
                        radius = 'NAN'
                        bins,fraction_cellType_0 = [0],[0]
                    
                    radius_table.write(growth_difference[gd]+','+density[de]+','+str(rep)+','+str(radius)+'\n')
                    for i in range(len(bins)):
                        radial_fraction_table.write(growth_difference[gd]+','+density[de]+','+str(rep)+','+str(bins[i])+','+str(fraction_cellType_0[i])+'\n')

                    
    radial_fraction_table.close()
    radius_table.close()
main()
