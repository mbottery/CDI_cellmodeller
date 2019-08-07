# import cellmodeller 
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
from CellModeller.gbStack import gbStack
import numpy
import math
import sys
import time
import random

'''
0-2-1
'''
# set up rates
inhibited_growth = 0.8
ratio = 0.5
density = 200
inhibitionRate = 0.1

type_0_cells = int(density*ratio)
type_3_cells = int(density*(1-ratio))

# set up sim
mu = 3.5
sigma = 0.5
max_cells = 200000
recoveryRate = 0.1

# cell types:
# 0 = CDI+
# 1 = CDI- in contact with CDI+
# 2 = CDI- inhibited
# 3 = CDI- not in contact with CDI+

growthRates = {0 : 1,
			   1 : 1,
			   2 : inhibited_growth,
               3 : 1}
dt = 0.05
seed = time.time()
random.seed(seed)
heap = gbStack(0,seed)
print("Seed: "+str(seed))

# set up colours
colours = {0 : [51/255.,153/255.,255/255.],
		   1 : [255/255.,102/255.,102/255.],
		   2 : [1,0,0],
           3 : [255/255.,153/255.,51/255.]}

def setup(sim):
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False,compNeighbours=True)
    regul = ModuleRegulator(sim, __file__)
    
    sim.init(biophys, regul, None, None)
    heap.clean()
    
    # randomly in a circle
    for i in xrange(type_0_cells):
        r = random.uniform(0.0, math.pi)
        x,y = circleRandPoint(0, 0, 200)
        sim.addCell(cellType=0, length=random.uniform(2, 5.8), pos=(x, y, 0.5), dir=(math.sin(r), math.cos(r), 0))
    for i in xrange(type_3_cells):
        r = random.uniform(0.0, math.pi)
        x,y = circleRandPoint(0, 0, 200)
        sim.addCell(cellType=3, length=random.uniform(2, 5.8), pos=(x, y, 0.5), dir=(math.sin(r), math.cos(r), 0))
    
    heap.addElement(dt,"step")

    if sim.is_gui:
        therenderer =Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)
    else:
        print "Running in batch mode: no display will be output"
    sim.savePickle = True
    sim.pickleSteps = 100


def init(cell):
    cell.targetVol = random.gauss(mu, sigma)
    cell.growthRate = growthRates[cell.cellType]
    cell.color = colours[cell.cellType]
    # propensities
    cell.reactions = {"inhibit":{"old":0,"current":0},
                      "recover":{"old":0,"current":0}}
    cell.divisionsInContact = 0
    cell.stepsInContact = 0


def update(cells,sim):
    # loop through all the cells in the simulation and calculate current
    # reaction rates
    for (id, cell) in cells.iteritems():
        update_cell(cell,cells)
        
    # given the current rates
    rxnID = heap.nextReaction().split('.')
    
    while rxnID[0] != "step":
        if rxnID[0] == "inhibit":
            # reaction is inhibit
            # need to find out which cell has become inhibited and change its
            # cell type to cell 2 and up date its rates in the model and on 
            # the stack
            id = sim.idxToId[int(rxnID[1])]
            # set its old and current rates to zero, choosing the reaction
            # removes it from the stack so if old > current it will produce
            # error that the id isn't on the stack as it tries to remove
            # it again
            cells[id].reactions['inhibit']['old'] = 0
            cells[id].reactions['inhibit']['current'] = 0
            cells[id].cellType = 2
            update_cell(cells[id],cells)
            # then we need to pick the next reaction on the stack and fire
            # that one
            rxnID = heap.nextReaction().split('.')
        if rxnID[0] == "recover":
            # reaction is recover 
            # need to find which cell has recovered from CDI and change its
            # cell type to 1 and update its rates in the model and on the stack
            id = sim.idxToId[int(rxnID[1])]
            # set its old and current rates to zero, choosing the reaction
            # removes it from the stack so if old > current it will produce
            # error that the id isn't on the stack as it tries to remove
            # it again
            cells[id].reactions['recover']['old'] = 0
            cells[id].reactions['recover']['current'] = 0
            cells[id].cellType = 3
            update_cell(cells[id],cells)
            # then we need to pick the next reaction on the stack and fire
            # that one
            rxnID = heap.nextReaction().split('.')
    heap.addElement(dt,"step");



def update_cell(cell,cells):
    # set divide flag if cell is long enough
    if cell.volume > cell.targetVol:
        cell.divideFlag = True
    
    if cell.cellType == 3:
        # if a CDI- cell comes into contact with a CDI+ cell set its type
        # to cellType 1 and set the division and step counter to zero
        if len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0]) > 0:
            cell.cellType = 1
            cell.divisionsInContact = 0
            cell.stepsInContact = 0
            
    # if the cell is non inhibited CDI-  
    if cell.cellType == 1:
        # if current cell is a CDI- cell in contact with a CDI+ cells, 
        # look at its contacts and check how many CDI+ cells it is 
        # touching and calculate inhibition rate based on the number of 
        # CDI+ contacts

        # linear increase in inhibition rate with increasing 
        # numbers of contacts, where 0 contacts = 0 inhibition rate

        # this loops through each of the cells contacts and creates
        # a list of all cells of type 0, the lenght of this list is the
        # number of CDI contacts the CDI- cell currently has
        CDI_contacts = len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0])
        
        cell.reactions['inhibit']['old'] = cell.reactions['inhibit']['current'] 
        cell.reactions['inhibit']['current'] = inhibitionRate*CDI_contacts
        # if the two rates are the same the cell has not gained or lost
        # and CDI+ contacts and the raction rate does not need to change
        # on the stack, otherwise, if the rate is different the reaction
        # needs adding, updating or removing from the stack
        
        if CDI_contacts == 0:
            # if the cell has lost contact with a CDI+ cell before 
            # becoming inhibited set its type back to cellType 3
            cell.cellType = 3
        else:
            # otherwise increase the step counter by 1
            cell.stepsInContact += 1
    if cell.cellType == 2:
        # Check to see if the inhibited cell is still touching a CDI+
        # cell, it is not set its recovery rate, this will add recovery
        # to the stack if the previous recovery rate was 0 and do nothing
        # if it previously had a recovery rate.
        if len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0]) == 0: 
            cell.reactions['recover']['old'] = cell.reactions['recover']['current']
            cell.reactions['recover']['current'] = recoveryRate
        # if the cell is in contact it sets the recovery rate to zero
        # if it used to have a recovery rate the reaction will be removed
        # from the stack, if it was previously in contact with a CDI+
        # cell this will do nothing.
        else:
            cell.reactions['recover']['old'] = cell.reactions['recover']['current']
            cell.reactions['recover']['current'] = 0
        
        
    # update colours and growthrates of the cell
    cell.growthRate = growthRates[cell.cellType]
    cell.color = colours[cell.cellType]

    # use the up-to-date rates to up-date the rates on the stack
    heapUpdateCell(cell)

def divide(parent, d1, d2):
    if parent.cellType == 1:
        d1.divisionsInContact += 1
        d2.divisionsInContact += 1
    d1.targetVol = random.gauss(mu, sigma)
    d2.targetVol = random.gauss(mu, sigma)
    
    # d1 already on stack, d2 set rates to zero, they will be calculated
    # next update and added to the stack if necessary
    d2.reactions["inhibit"]["old"] = 0
    d2.reactions["inhibit"]["current"] = 0
    d2.reactions["recover"]["old"] = 0
    d2.reactions["recover"]["current"] = 0   
    

def heapUpdateCell(cell):
    # sort out rates on stack:
    # if old rate is 0 and current rate > 0 need to add it to the stack
    # if old rate is > 0  and current rate = 0 need to remove the reaction
    # else if old rate != current rate need to update the reaction 
    for (id, r) in cell.reactions.iteritems():
        rxnID = id+'.'+str(cell.idx)
        if r['old'] == 0 and r['current'] > 0:
            heap.addElement(-1.0*math.log(random.random())/r['current'],rxnID, r['current'])
        elif r['old']  > 0 and r['current'] == 0:
            heap.removeOne(rxnID)
        elif r['old']  != r['current']:
            heap.updateElement(rxnID,r['old'] , r['current'])

	
def circleRandPoint(x1,y1,rc):
    a = 2*math.pi*random.random()
    r = math.sqrt(random.random())
    x = (rc*r)*math.cos(a)+x1
    y = (rc*r)*math.sin(a)+y1
    return(x, y)


