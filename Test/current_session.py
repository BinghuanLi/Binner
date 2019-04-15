# coding: utf-8
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import KDTree

Points = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
Vor = Voronoi(Points)
Vor.vertices
Vor.points
Vor.point_region
Vor.regions
tree = KDTree(Points)
tree.query([1, 0.1])

#get_ipython().magic(u'save current_session ~0/')
