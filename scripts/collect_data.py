'''
Counts the number of cells through time
1) in the total population
2) of each cell type
'''
import sys
import os
sys.path.append('.')
import cPickle
import glob

def openPickle(pickle):
    print('opening: '+ pickle)
    data = False
    try:
        data = cPickle.load(open(pickle,'r'))
    except:
        print 'Failed to open pickle file'
        print sys.exc_info()[0]
    return data
        
def countCells(data):
    cs = data['cellStates']
    it = iter(cs)
    n = len(cs)
    n0 = 0
    n1 = 0
    n2 = 0
    n3 = 0
    for it in cs:
        if cs[it].cellType == 0:
            n0+=1
        if cs[it].cellType == 1:
            n1+=1
        if cs[it].cellType == 2:
            n2+=1
        if cs[it].cellType == 3:
            n3+=1
    return n,n0,n1,n2,n3

def main():
    # open all of the step files in each of the directories
    # count number of cells of each cell type
    # if pickle step == 700 count the number of cells in the centre of the colony
    # return this data to a main data table with all the other pickle files
    # make sure the table has all the parameters and pc which ran the simulation 
    # and replicate number
    # write out the table so I can plot it in R
    indirs = sys.argv[1:]
    
    countfile = open('growth_difference_count.txt','w')
    countfile.write('growth_difference,density,rep,step,step*dt,n,n0,n1,n2,n3\n')
    
    growth_difference = {'0':'0','1':'0.01','2':'0.05','3':'0.1','4':'0.2'}
    density = {'0':'20000','1':'2000','2':'200'}
    
    codes = []
    
    for d in range(len(indirs)): # machines
        for dirs in os.listdir(indirs[d]):  # simulations 
            print os.path.join(indirs[d],dirs)
            # get the parameter set from the directory name
            # add it to the code list and figure out rep 
            # number
            code = dirs[18:21]
            print 'simulation: {0}'.format(code)
            codes.append(code)
            rep = indirs[d][3:]
            gd,de = code[0],code[2]
            print 'rep: {0}'.format(rep)
            for filename in os.listdir(os.path.join(indirs[d],dirs)): # pickle steps
                if filename == "step-00700.pickle":
                    step = int(filename.split('.')[0].split('-')[1])
                    data = openPickle(os.path.join(indirs[d],dirs,filename))
                    if data:
                        n,n0,n1,n2,n3 = countCells(data)
                    else:
                        n,n0,n1,n2,n3 = ['NAN','NAN','NAN','NAN','NAN']
                    countfile.write(growth_difference[gd]+','+density[de]+','+str(rep)+','+str(step)+','+str(step*0.05)+','+str(n)+','+str(n0)+','+str(n1)+','+str(n2)+','+str(n3)+'\n')
    
    countfile.close()
main()

