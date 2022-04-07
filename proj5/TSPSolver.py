#!/usr/bin/python3

import copy
from re import L
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		current_city_idx = 0
		start_time = time.time()
		while current_city_idx < ncities and not foundTour and time.time()-start_time < time_allowance: # O(ncities), worst case scenario, final route should start with the last city
			count += 1
			route = [cities[current_city_idx]]
			not_visited = set(range(ncities))
			not_visited.remove(current_city_idx)

			# O(ncities)
			# - when a solution is possible, the route should contain ncities
			# - when break is called, in worst case scenario, we can construct a route of ncities - 1, and then find out it's not possible to go to the remaining city
			while len(route) != ncities and time.time()-start_time < time_allowance:
				# from last city in route
				# find the city with minimum cost to
				min_cost = float('inf')
				next_city_idx = -1
				# O(ncities)
				# not_visited has n-1 cities at 1st iteration of the while loop, then n-2 at 2nd, n-3 ... 1 at last
				for city_idx in not_visited:
					cost = route[-1].costTo(cities[city_idx])
					if cost < min_cost:
						min_cost = cost
						next_city_idx = city_idx

				if next_city_idx == -1:
					break
				else:
					not_visited.remove(next_city_idx)
					route.append(cities[next_city_idx])

			if len(route) == ncities:
				if route[-1].costTo(route[0]) < float('inf'):
					bssf = TSPSolution(route)
					foundTour = True
	
			if not foundTour:
				current_city_idx += 1
		# final time complexity O(ncities ^ 3)
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None

		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		pass



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		bssf = None
		start_time = time.time()

		subsets = self.genPowerset(ncities)
		routeLenMap = [[math.inf for i in range(ncities)] for j in range (pow(2,ncities-1))]
		routeLenMap[0][0] = 0
		routeMap = [[[0] for i in range(ncities)] for j in range (pow(2,ncities-1))]
		for i in range(pow(2,ncities-1)):
			if (i == 0):
				routeLenMap[0][i] = 0
			else:
				routeLenMap[i][0] = math.inf
		for i in range(len(subsets)): #for all subsets
			self.getMin(routeMap, routeLenMap, subsets, subsets[i], i, cities)
		shortestPath = self.findSmallest(routeLenMap,pow(2,ncities-1)-1, cities, 0)
		#return values for results
		end_time = time.time()
		route = self.generateRoute(routeMap[pow(2,ncities-1)-1][shortestPath], cities)
		bssf = TSPSolution(route)
		results['cost'] = bssf.cost
		results['time']  = end_time - start_time
		results['count'] = 1
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results
		#TODO: time

	def generateRoute(self, currRoute, cities):
		route = []
		#for i in range(len(currRoute)-1, -1, -1):
		for i in range(len(currRoute)):
			route.append(cities[currRoute[i]])
		return route

	def genPowerset(self,ncities):
		cities = [i for i in range(ncities)]
		subsets = [0 for i in range(pow(2,ncities-1))]
		index = 1
		firstTuple = (0,0)
		subsets[0] = firstTuple
		for i in range(2, len(cities)):
			for element in itertools.combinations(cities,i):
				if (element.count(0) == 1):
					subsets[index] = element
					index += 1
		lastTuple = tuple(cities)
		subsets[len(subsets)-1] = lastTuple
		return subsets
	
	def getMin(self, routeMap, routeLenMap, subsets, subset, index, cities): #need to change what i we are getting, 
		#currently it is which subset we are on and we want it to be where we are trying to endup at
		dist = math.inf
		for j in range(len(subset)): #for all j in subset where j !=0
			if (subset[j] != 0):
				if (len(subset) > 2):
					oldDistSubset = subset[:j] + subset[j+1:] #gets the subset - j
					subsetIndex = subsets.index(oldDistSubset) #gets what index that would be on our table
					lastCity = self.findSmallest(routeLenMap,subsetIndex, cities, subset[j]) #finds which last city had the shortest path
					oldDist = routeLenMap[subsetIndex][lastCity] #gets the dist of the last route we took
					newRoute = copy.deepcopy(routeMap[subsetIndex][lastCity]) #gets the array of cities on our route
				else:
					oldDist = 0
					lastCity = 0
					newRoute = copy.deepcopy(routeMap[index][0])
				dist = oldDist + cities[lastCity].costTo(cities[subset[j]]) #distance from last city on route to next city
				newRoute.append(subset[j])
				routeMap[index][subset[j]] = newRoute
				routeLenMap[index][subset[j]] = dist			

	def findSmallest(self, routeLenMap, subsetIndex, cities,goalCity):
		value = math.inf
		index = 0
		for i in range(len(cities)):
			oldDist = routeLenMap[subsetIndex][i]
			nextDist = cities[i].costTo(cities[goalCity])
			if (oldDist+nextDist < value):
				value = oldDist + nextDist
				index = i
		return index

	