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
1-1-0-0
'''
# set up rates
inhibited_growth = 0.0
density = 2000
inhibitionRate = 1
cost = 0

ratio = 0.5

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

growthRates = {0 : 1-(1*cost),
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
    for (id, cell) in cells.iteritems():
        update_cell(cell,cells)
    
    rxnID = heap.nextReaction().split('.')
    
    while rxnID[0] != "step":
        if rxnID[0] == "inhibit":
            id = sim.idxToId[int(rxnID[1])]

            cells[id].reactions['inhibit']['old'] = 0
            cells[id].reactions['inhibit']['current'] = 0
            cells[id].cellType = 2
            update_cell(cells[id],cells)

            rxnID = heap.nextReaction().split('.')
        if rxnID[0] == "recover":

            id = sim.idxToId[int(rxnID[1])]

            cells[id].reactions['recover']['old'] = 0
            cells[id].reactions['recover']['current'] = 0
            cells[id].cellType = 3
            update_cell(cells[id],cells)

            rxnID = heap.nextReaction().split('.')

    heap.addElement(dt,"step");



def update_cell(cell,cells):
    if cell.volume > cell.targetVol:
        cell.divideFlag = True
    
    if cell.cellType == 3:
        if len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0]) > 0:
            cell.cellType = 1
            cell.divisionsInContact = 0
            cell.stepsInContact = 0
              
    if cell.cellType == 1:
        CDI_contacts = len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0])
        
        cell.reactions['inhibit']['old'] = cell.reactions['inhibit']['current'] 
        cell.reactions['inhibit']['current'] = inhibitionRate*CDI_contacts

        if CDI_contacts == 0:
            cell.cellType = 3
        else:
            cell.stepsInContact += 1
    if cell.cellType == 2:
        if len([cells[index].cellType for index in cell.neighbours if cells[index].cellType == 0]) == 0: 
            cell.reactions['recover']['old'] = cell.reactions['recover']['current']
            cell.reactions['recover']['current'] = recoveryRate
        else:
            cell.reactions['recover']['old'] = cell.reactions['recover']['current']
            cell.reactions['recover']['current'] = 0

    cell.growthRate = growthRates[cell.cellType]
    cell.color = colours[cell.cellType]

    heapUpdateCell(cell)

def divide(parent, d1, d2):
    if parent.cellType == 1:
        d1.divisionsInContact += 1
        d2.divisionsInContact += 1
    d1.targetVol = random.gauss(mu, sigma)
    d2.targetVol = random.gauss(mu, sigma)
    
    d2.reactions["inhibit"]["old"] = 0
    d2.reactions["inhibit"]["current"] = 0
    d2.reactions["recover"]["old"] = 0
    d2.reactions["recover"]["current"] = 0   
    

def heapUpdateCell(cell):
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


