import numpy as np
import sys
import os
sys.path.append('.')
import cPickle
import glob
import operator
import matplotlib.pyplot as plt

# split the cells up into annuli

# collaps the cells down on a line in the centre of the annuli

# equation of a circle: x^2 + y^2 = r^2

'''
where P is the point, C is the center, and R is the radius
double vX = pX - cX;
double vY = pY - cY;
double magV = sqrt(vX*vX + vY*vY);
double aX = cX + vX / magV * R;
double aY = cY + vY / magV * R;
'''


class radial_fft():
	# Class that collects radial data from pickle files
	def __init__(self,data,bins,p=False):
		self.states = data['cellStates']
		self.bins = bins

		# make the empty data frames
		self.data = []
		self.radialPos = []
		self.annuli=[[] for i in range(len(bins)-1)]
		self.collapsed_annuli=[[] for i in range(len(bins)-1)]
		self.dist=[[] for i in range(len(bins)-1)]
		self.angle=[[] for i in range(len(bins)-1)]
		self.arc_dist=[[] for i in range(len(bins)-1)]
		# bin the data
		self.bin_data_annuli()
		# sort each circle of points 
		for i in range(len(self.collapsed_annuli)):
			self.sort_circle_points(self.collapsed_annuli[i])
		# calculate angle and distas
		self.angle_distance()
		
		if p:
			self.plot_data()

	def calc_radius(self):
		return max(self.radialPos)
    	
	def bin_data_annuli(self):
		# loop through each cell
		for (id,s) in self.states.iteritems():
			radial_position = np.sqrt(s.pos[0]**2+s.pos[1]**2)
			self.radialPos.append(radial_position)
			self.data.append([s.pos[0],s.pos[1],s.cellType])
			# loop through each bin
			for j in range(len(self.bins)-1):	
				# if radius of point 
				if self.bins[j]<=radial_position<self.bins[j+1]:
					# add cells what lie within the annulus to the appropriate bin
					self.annuli[j].append(self.data[-1])
					# calculate the centre of the annulus
					R = (self.bins[j]+self.bins[j+1])/2.
					# project points onto the centre line
					projected_point=list(self.project_points_to_radius(R,self.data[-1][0],self.data[-1][1]))
					# add on the cell type to new coordinate
					projected_point.append(self.data[-1][2])
					# add this to the data
					self.collapsed_annuli[j].append(projected_point)

	def angle_distance(self):
		for i in range(len(self.collapsed_annuli)):
			# calculate the current radius of the point (the mid point of the current annulus)
			R = (self.bins[i]+self.bins[i+1])/2.
			# loop through each point
			for j in range(len(self.collapsed_annuli[i])):
				# if it isn't the last point
				if j != len(self.collapsed_annuli[i])-1:
					# calculate the distance between this point and the next
					self.dist[i].append(self.distance_between_points(self.collapsed_annuli[i][j][0],self.collapsed_annuli[i][j][1],self.collapsed_annuli[i][j+1][0],self.collapsed_annuli[i][j+1][1]))
					# calculate the angle between this point and the next using the distance
					self.angle[i].append(self.angle_between_two_points(self.dist[i][-1],R))
					# calculate the arc distance between the two points
					self.arc_dist[i].append(self.angle[i][-1]*R)
		
				# if it is the last point 
				else:
					# calculate the the distance between this point and the first point
					self.dist[i].append(self.distance_between_points(self.collapsed_annuli[i][j][0],self.collapsed_annuli[i][j][1],self.collapsed_annuli[i][0][0],self.collapsed_annuli[i][0][1]))
					# calculate the angle between this point and the first using the distance
					self.angle[i].append(self.angle_between_two_points(self.dist[i][-1],R))
					# calculate the arc distance between the two points
					self.arc_dist[i].append(self.angle[i][-1]*R)

	def distance_between_points(self,x1,y1,x2,y2):
		# finds the distance between two points in 2D space
		d = np.sqrt((x1-x2)**2+(y1-y2)**2)
		return d

	def angle_between_two_points(self,d,r):
		# calculates the angle in radians between two points on a circle with the radius R
		angle = np.arccos(1-((d**2)/(2*(r**2))))
		return angle
		
	def project_points_to_radius(self,R,X,Y):
		# takes a point in 2D space and finds closest point on circle with radius R
		# with origin 0,0
		magV = np.sqrt(X**2 + Y**2)
		aX = X/magV*R
		aY = Y/magV*R
		return aX,aY

	def sort_circle_points(self,c):
		# order the points from a given position, in this case atan does it from the x 
		# axis clockwise
		c.sort(key=lambda c:np.arctan2(c[0], c[1]))

	def plot_data(self):
		for i in range(len(self.annuli)):
			R = (self.bins[i]+self.bins[i+1])/2.
			fig, ax = plt.subplots(1,3) 
			ax[0].scatter([item[0] for item in self.data],[item[1] for item in self.data],c=[item[2] for item in self.data])
			ax[0].set_xlim((-400,400))
			ax[0].set_ylim((-400,400))
			circle = plt.Circle((0, 0), R, color = 'r', fill=False)
			ax[0].add_artist(circle)
			for j in range(len(self.bins)):
				circle1 = plt.Circle((0, 0), self.bins[j], color = 'b', fill=False)
				ax[0].add_artist(circle1)
			
			ax[1].scatter([item[0] for item in self.annuli[i]],[item[1] for item in self.annuli[i]],c=[item[2] for item in self.annuli[i]])
			ax[1].set_xlim((-400,400))
			ax[1].set_ylim((-400,400))
			circle = plt.Circle((0, 0), R, color = 'r', fill=False)
			ax[1].add_artist(circle)
			for j in range(len(self.bins)):
				circle1 = plt.Circle((0, 0), self.bins[j], color = 'b', fill=False)
				ax[1].add_artist(circle1)
			
			ax[2].scatter([item[0] for item in self.collapsed_annuli[i]],[item[1] for item in self.collapsed_annuli[i]],c=[item[2] for item in self.collapsed_annuli[i]])
			ax[2].set_xlim((-400,400))
			ax[2].set_ylim((-400,400))
			circle = plt.Circle((0, 0), R, color = 'r', fill=False)
			ax[2].add_artist(circle)
			for j in range(len(self.bins)):
				circle1 = plt.Circle((0, 0), self.bins[j], color = 'b', fill=False)
				ax[2].add_artist(circle1)
			plt.show()
	
	def output_data(self):
        	radii = [(self.bins[i]+self.bins[i+1])/2. for i in range(len(self.bins)-1)]
		return radii,self.collapsed_annuli,self.arc_dist
    
	def write_data(self,fileName):
		df = open(fileName,'w')
		df.write('radius,cellType,distance\n')
		for i in range(len(self.annuli)):
			R = (self.bins[i]+self.bins[i+1])/2.
			for j in range(len(self.arc_dist[i])):
				df.write(str(R)+','+str(self.collapsed_annuli[i][j][2])+','+str(self.arc_dist[i][j])+'\n')
	      
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
    indirs = sys.argv[1:]
    
    
    radial_fft_table = open('growth_radial_fft.txt','w')
    radial_fft_table.write('growth_inhibited,density,rep,radius,cellType,distance\n')
    
    growth_difference = {'0':'0','1':'0.01','2':'0.05','3':'0.1','4':'0.2'}
    density = {'0':'20000','1':'2000','2':'200'}
    
    codes=[]
    
    for d in range(len(indirs)): # machine
    	print indirs[d]
        for dirs in os.listdir(indirs[d]):  # simulations 
            print os.path.join(indirs[d],dirs)
            rep = indirs[d][3:]
            code = dirs[18:21]
            print 'simulation: {0}'.format(code)
            codes.append(code)
            gd,de = code[0],code[2]
            print 'rep: {0}'.format(rep)
            for filename in os.listdir(os.path.join(indirs[d],dirs)): # pickle steps
                if filename == "step-00700.pickle":
                    data = openPickle(os.path.join(indirs[d],dirs,filename))
                    if data:
                        r = radial_fft(data,range(200,420,20),p=False)
                        radii,cell_types,distances = r.output_data()
                    else:
                        # blank data
                        print('no data')
                        radii,cell_types,distances = [0],[0],[0]
                    for i in range(len(radii)):
                        for j in range(len(distances[i])):
                            radial_fft_table.write(growth_difference[gd]+','+density[de]+','+str(rep)+','+str(radii[i])+','+str(cell_types[i][j][2])+','+str(distances[i][j])+'\n')
			
      
    radial_fft_table.close()
    
main()

