#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import sys 

class Node:

	def __init__(self, ID, vertices):

		self.ID = ID 
		self.vertices = vertices 

		self.boundary = 0 

	def setBoundary(self, val=1):

		self.boundary = val 

	def getString(self):

		string = '{} {} {} {} {}\n'.format(self.ID-1, self.boundary, 
			self.vertices[0], self.vertices[1], self.vertices[2])

		return string 

class Element:

	def __init__(self, ID, node_list, group=0):

		self.ID = ID 
		self.node_list = node_list
		self.group = group

	def getString(self, skip, gskip):

		string = str(self.ID-skip-1) + ' ' + str(self.group-gskip) + ' '
		for i in range(len(self.node_list)):
			string += str(self.node_list[i]-1) + ' '

		return string + '\n'

if (len(sys.argv) < 2):
	print('specify file name')
	sys.exit()

fname = sys.argv[1]

f = open(fname, 'r')

mnodes = [] 
mele = [] 

boundary_elements = 0 

for line in f:

	if (line.startswith('$Nodes')):

		nnodes = int(f.readline().strip())

		for i in range(nnodes):

			line = f.readline().strip().split(' ')

			ID = int(line[0])
			vertices = line[1:]
			for j in range(3):
				vertices[j] = float(vertices[j])

			mnodes.append(Node(ID, vertices))

	elif (line.startswith('$Elements')):

		nel = int(f.readline().strip())

		for i in range(nel):

			line = f.readline().strip().split(' ')

			ID = int(line[0]) 

			# first order boundary line 
			if (int(line[1]) == 1):

				skip = int(line[2]) + 3
				nn = int(line[skip])
				mnodes[nn-1].setBoundary(int(line[3]))
				nn = int(line[skip+1])
				mnodes[nn-1].setBoundary(int(line[3]))

				boundary_elements += 1

			# second order boundary line
			elif (int(line[1]) == 8):
				nodes = line[len(line)-3:]
				for j in range(len(nodes)):
					val = int(nodes[j])
					mnodes[val-1].setBoundary(int(line[3]))
				boundary_elements += 1

			# third order boundary line 
			elif (int(line[1]) == 26):
				nodes = line[len(line)-4:]
				for j in range(4):
					val = int(nodes[j])
					mnodes[val-1].setBoundary(int(line[3]))
				boundary_elements += 1

			# first order element 
			elif (int(line[1]) == 2):
				nodes = line[len(line)-3:]
				for j in range(len(nodes)):
					nodes[j] = int(nodes[j])
				mele.append(Element(ID, nodes, int(line[4])))

			# second order element 
			elif (int(line[1]) == 9):
				nodes = line[len(line)-6:]
				for j in range(6):
					nodes[j] = int(nodes[j])
				mele.append(Element(ID, nodes, int(line[4])))

			# third order element 
			elif (int(line[1]) == 21):
				nodes = line[len(line)-10:]
				for j in range(10):
					nodes[j] = int(nodes[j])
				mele.append(Element(ID, nodes, int(line[4])))

			else:
				print('element type not defined')
				sys.exit()

f.close()

minGroup = 10000
for i in range(len(mele)):
	if (mele[i].group < minGroup):
		minGroup = mele[i].group 

nodefile = open('.'.join(fname.split('.')[:-1])+'.node', 'w')
nodefile.write(str(len(mnodes)) + '\n')
for i in range(len(mnodes)):

	nodefile.write(mnodes[i].getString())

nodefile.close()
			
elfile = open('.'.join(fname.split('.')[:-1])+'.ele', 'w')
elfile.write(str(len(mele)) + '\n')
for i in range(len(mele)):

	elfile.write(mele[i].getString(boundary_elements, minGroup))

elfile.close()

# NODES = mele[0].node_list
# for n in NODES:

# 	xy = mnodes[n-1].vertices

# 	print(xy)

# 	plt.plot(xy[0], xy[1], 'o')

# plt.show()