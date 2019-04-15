# coding: utf-8
import numpy as np
import pandas as pd
import tool
from tool import load_data_2017
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

# set mimimum bkg for significance calculation
minimum_bkg = 0.000001

# load data
inputPath = "../data/"
variables = ["mvaOutput_2lss_ttV","mvaOutput_2lss_ttbar"]
data=load_data_2017(inputPath, variables, "passGenMatchCut==1")

# load signal points to generate voronoi diagram
sig_points_x = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttV"]
sig_points_y = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttbar"]

Points=np.vstack((sig_points_x,sig_points_y)).T

Vor = Voronoi(Points)
print("print voronoi vertices")
print(Vor.vertices)
print("print voronoi regions")
print(Vor.regions)

# set signal vor_point and vor_region
mask = data.target ==1
data.loc[mask,"vor_region"]=Vor.point_region
data.loc[mask,"vor_point"]=np.arange(len(Vor.points))

#fig = voronoi_plot_2d(Vor, show_vertices=False, line_colors='orange', line_width=2, line_alpha=0.6, point_size=2)

# set ttV vor_point and vor_region
ttV_x = data.ix[(data.target.values==2)]["mvaOutput_2lss_ttV"]
ttV_y = data.ix[(data.target.values==2)]["mvaOutput_2lss_ttbar"]
ttV_points = np.vstack((ttV_x,ttV_y)).T

ttJ_x = data.ix[(data.target.values==3)]["mvaOutput_2lss_ttV"]
ttJ_y = data.ix[(data.target.values==3)]["mvaOutput_2lss_ttbar"]
ttJ_points = np.vstack((ttJ_x,ttJ_y)).T

point_idx, region_idx = tool.find_vor_region(ttV_points,Vor)
print("region idx is {}".format(region_idx))
mask = data.target ==2
data.loc[mask,"vor_region"]= region_idx
data.loc[mask,"vor_point"]= point_idx

point_idx, region_idx = tool.find_vor_region(ttJ_points,Vor)
print("region idx is {}".format(region_idx))
mask = data.target ==3
data.loc[mask,"vor_region"]= region_idx
data.loc[mask,"vor_point"]= point_idx

# sort data by vor_region
print (data)

# plot initial Voronoi diagrams
plt = tool.plot_vor_2d(Points, Vor.vertices, Vor.regions, Vor.ridge_vertices, Vor.ridge_points) 
plt.xlabel("mvaOutput_2lss_ttV")
plt.ylabel("mvaOutput_2lss_ttbar")
plt.plot(ttV_x,ttV_y,'r.', label='ttV')
plt.plot(ttJ_x,ttJ_y,'k.', label='ttJ')
plt.legend()
#plt.show()
plt.savefig("Voronoi_test.png")

# create df to save the initial voronoi cell informations
column_list = ["vor_point","vor_region","s","b1","b2","error","significance"]
df_vor =  pd.DataFrame(columns=column_list)
for cell in range(1+(data["vor_region"].max())):
    if cell not in data["vor_region"].values: continue
    v_point = data[(data["target"]==1) & (data["vor_region"]==cell)]["vor_point"].values[0]
    #print(v_point)
    v_region = cell
    nS = data[(data["target"]==1) & (data["vor_region"]==cell)]["totalWeight"].sum() 
    nB1 = data[(data["target"]==2) & (data["vor_region"]==cell)]["totalWeight"].sum() 
    nB2 = data[(data["target"]==3) & (data["vor_region"]==cell)]["totalWeight"].sum() 
    err_b1 = data.ix[data.target.values==2]["error"].values[0]
    err_b2 = data.ix[data.target.values==3]["error"].values[0]
    err_sum_b , significance = tool.get_significance(nS, nB1, nB2, err_b1, err_b2, minimum_bkg)
    a = [v_point,  v_region, nS, nB1,nB2,err_sum_b,significance]
    df = pd.DataFrame([a],columns=column_list)
    df_vor = df_vor.append(df, ignore_index=True)

print(df_vor)
