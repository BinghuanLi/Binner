# coding: utf-8
import numpy as np
import pandas as pd
import tool
import math
from tool import load_data_2017
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from colorized_voronoi import voronoi_finite_polygons_2d
from shapely import geometry
import shapely
import geopandas


# set mimimum bkg for significance calculation
minimum_bkg = 0.000001
negative_weight = "Sort" # Sort/Dist
delta = 1000 # load event every delta instance, so the events loaded would be N/delta
group_events = -1 # group events so each initial cell contains group_events instance, the initial #cell = (N/delta)/group_events
doPlot = True
load_csv = True
err_b1 = 0.15 
err_b2 = 0.25 

# load data
inputPath = "../data/"
variables = ["mvaOutput_2lss_ttV","mvaOutput_2lss_ttbar"]
if not load_csv:
    data=load_data_2017(inputPath, variables, "passGenMatchCut==1", skip = delta ,nEvt=group_events, pair_negwgt = negative_weight ) 
    data.to_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}.csv".format(inputPath,negative_weight,delta,group_events))
else:
    try: data = pd.read_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}.csv".format(inputPath,negative_weight,delta,group_events))
    except:
        print( "{}input_PairNegWgt{}_Delta{}_GroupEvt{}.csv deosn't exist ".format(inputPath,negative_weight,delta,group_events))
        print( "run load_data_2017 ")
        data=load_data_2017(inputPath, variables, "passGenMatchCut==1", skip = delta ,nEvt=group_events, pair_negwgt = negative_weight ) 
        data.to_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}.csv".format(inputPath,negative_weight,delta,group_events))

print (data)

# save log file
#file_log = open(("Significance_Initial_PairNegWgt{}_Delta{}_GroupEvt{}.log".format(negative_weight,delta,group_events)),"w")
#file_log.write(df_vor)
#file_log.close()

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
print ("finish assigning bkgs")
print (data)

if doPlot and len(Vor.regions) < 1000 :
#if doPlot:
    # plot initial Voronoi diagrams only if vor.regions < 1000
    print (" plot initial Voronoi diagrams ")
    plt = tool.plot_vor_2d(Points, Vor.vertices, Vor.regions, Vor.ridge_vertices, Vor.ridge_points) 
    #plt.xlabel("mvaOutput_2lss_ttV")
    #plt.ylabel("mvaOutput_2lss_ttbar")
    #plt.plot(ttV_x,ttV_y,'r.', label='ttV')
    #plt.plot(ttJ_x,ttJ_y,'k.', label='ttJ')
    plt.legend()
    plt.title("initial voronoi diagrams")
    #plt.show()
    plt.savefig("Voronoi_initial_PairNegWgt{}_Delta{}_GroupEvt{}.png".format(negative_weight,delta,group_events))
    plt.close()

    print (" finish plot initial Voronoi diagrams ")

#column_list = ["vor_point","vor_region","s","b1","b2","b1_err","b2_err","error","significance"]
my_values = ["vor_point","totalWeight","error"]
my_index = ["vor_region","target"]
#print(df)
my_list_column = my_values 
my_list_agg = ['mean']*len(my_list_column)
my_dict = dict(zip(my_list_column,my_list_agg))
my_dict['totalWeight']='sum'
data_table = data.pivot_table(values=my_values, index=my_index, aggfunc = my_dict )
data_table = data_table.unstack('target', fill_value=0).stack(level=1)
print data_table
my_vor_region = data_table.index.get_level_values(0)[data_table.index.get_level_values(1)==1].values
print (my_vor_region)
my_vor_point = data_table.loc[data_table.index.get_level_values(1)==1]["vor_point"].values
print (my_vor_point)
my_s = data_table.loc[data_table.index.get_level_values(1)==1]["totalWeight"].values
my_b1 = data_table.loc[data_table.index.get_level_values(1)==2]["totalWeight"].values
my_b2 = data_table.loc[data_table.index.get_level_values(1)==3]["totalWeight"].values
my_err_sum_b = np.zeros(len(my_vor_region))
my_significance = np.zeros(len(my_vor_region))
print(my_err_sum_b)
print ("calculate significance")
vfunc = np.vectorize(tool.get_significance)
my_err_sum_b, my_significance = vfunc(my_s, my_b1, my_b2, err_b1, err_b2, minimum_bkg)

print ("finish calculate significance")

df_vor =  pd.DataFrame({"vor_point":my_vor_point,"vor_region":my_vor_region,"s":my_s,"b1":my_b1,"b2":my_b2,"error":my_err_sum_b,"significance":my_significance})

print(df_vor)

# sort by significance
df_vor.sort_values(by="significance", ascending=False, inplace=True)
df_vor = df_vor.reset_index(drop=True)

print(df_vor)
print(df_vor["significance"].values)

if doPlot:
    # plot the original significance
    print (" plot original significance ")
    values1, bins, _ = plt.hist(df_vor["significance"].values, bins = 10, alpha =0.5, range=(0., 1.01*df_vor["significance"].values[0]) , density = True)
    plt.xlabel("Z")
    plt.ylabel("%")
    plt.title("Z distribution")
    plt.savefig("Significance_Initial_PairNegWgt{}_Delta{}_GroupEvt{}.png".format(negative_weight,delta,group_events))
    plt.close()

# bin merge algorithm
s_trail =0
b1_trail = 0
b2_trail = 0
Z_trail = 0
n_cells = len(df_vor.index)
# set cost function
cost_fom = 0
sum_Zsquare =  (df_vor["significance"]**2).sum()
sum_deltaZ =0
# save updated Z
total_Z = np.zeros(n_cells)
# save Z before updates
lag_Z = np.zeros(n_cells)
# save flag to indicate wether merge happens
update_Z = np.zeros(n_cells)
# save the new vor_cell number, starting from 0
vor_cell = np.full((n_cells,),-1)
vor_count = 0
for index, row in df_vor.iterrows():
    s = row["s"]
    b1 = row["b1"]
    b2 = row["b2"]
    Z = row["significance"]
    lag_Z_square = sum_Zsquare
    lag_delta_Z = sum_deltaZ
    err_sum_b , significance = tool.get_significance((s_trail+s), (b1_trail+b1), (b2_trail+b2), err_b1, err_b2, minimum_bkg)
    #print("s_trail:{}, b1_trail:{}, b2_trail:{}, s:{}, b1:{}, b2:{}, Cost:{}, Significance: {}, Sum_Zsqure: {}".format(s_trail, b1_trail, b2_trail, s, b1, b2, cost,(significance**2),(Z**2+Z_trail**2)))
    cost = tool.cost_fun(significance, Z_trail, Z, cost_fom)
    if cost >= 0.:
        # merge the bins
        sum_Zsquare += cost
        sum_deltaZ += cost
        s_trail += s
        b1_trail += b1_trail
        b2_trail += b2_trail
        Z_trail = significance
        update_Z[index] = 1
        if cost_fom ==0:
            lag_Z[index] = math.sqrt(lag_Z_square)
            total_Z[index] = math.sqrt(sum_Zsquare)
        else:
            lag_Z[index] = sum_deltaZ
            total_Z[index] = sum_deltaZ
    else:
        s_trail = s
        b1_trail = b1
        b2_trail = b2
        Z_trail = Z
        if cost_fom ==0:
            lag_Z[index] = math.sqrt(lag_Z_square)
            total_Z[index] = math.sqrt(sum_Zsquare)
        else:
            lag_Z[index] = sum_deltaZ
            total_Z[index] = sum_deltaZ
        vor_count +=1
    vor_cell[index] = vor_count 

df_vor["total_sig"] = total_Z
df_vor["lag_sig"] = lag_Z
df_vor["update"] = update_Z
df_vor["vor_label"] = vor_cell

df_vor.to_csv("PairNegWgt{}_Delta{}_GroupEvt{}.csv".format(negative_weight,delta,group_events))
print(df_vor)

#if doPlot and len(Vor.regions) < 1000 :
if doPlot :
    # plot merged Voronoi diagrams with color
    print (" plot final merged Voronoi diagrams ")
    
    regions, vertices = voronoi_finite_polygons_2d(Vor.points, Vor.vertices, Vor.regions, Vor.ridge_vertices, Vor.ridge_points, Vor.point_region)

    final_label_size = df_vor["vor_label"].nunique()
    final_labels=df_vor.vor_label.unique()
    label_points = df_vor.groupby('vor_label')['vor_point'].apply(list).to_dict() 
    '''
    colors = cm.rainbow(np.linspace(0,1,final_label_size))
    points_indices = np.arange(len(regions))
    for region, p in zip(regions, points_indices):
        #print p, region
        m_label = df_vor[df_vor["vor_point"]==p]["vor_label"].values[0]
        c = colors[m_label]
        polygon = vertices[region]
        #print polygon
        plt.fill(*zip(*polygon), alpha=0.4, color=c)
    '''
    cmap = cm.ScalarMappable(cmap='rainbow')
    colors = cmap.to_rgba(final_labels)
    fig, ax = plt.subplots()
    print (" ready to loop over labels ")
    for label, points in label_points.iteritems():
        #print ("label: {}".format(label))
        if label % 1000 ==0: print ("label: {}".format(label))
        c = colors[label]
        geoPolys = []
        n_polygons = len(points)
        for p in points:
            polygon = vertices[regions[int(p)]] 
            geoPoly =  geometry.Polygon(polygon)
            geoPolys.append(geoPoly)
        s = geopandas.GeoSeries(geoPolys) 
        col = tool.plot_polygon_collection(ax, s.geometry, facecolor=c)
    print (" finish plot loop over labels ")

    #plt.plot(Points[:,0], Points[:,1], '.', label='ttH')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.title("Final Voronoi Diagram")
    plt.xlabel("mvaOutput_2lss_ttV")
    plt.ylabel("mvaOutput_2lss_ttbar")
    #plt.plot(ttV_x,ttV_y,'r.', label='ttV')
    #plt.plot(ttJ_x,ttJ_y,'k.', label='ttJ')
    plt.legend()

    #plt.show()
    plt.savefig("ColoredVoronoi_Final_PairNegWgt{}_Delta{}_GroupEvt{}.png".format(negative_weight,delta,group_events))
    plt.close()

if doPlot: 
    # plot the evolving history
    print (" plot the evolving history ")
    my_slice = int(math.ceil(n_cells/100.))
    evolve_x = np.arange(n_cells)
    y_totsig = df_vor["total_sig"].values
    y_lagsig = df_vor["lag_sig"].values
    y_deltasig = 100.*(y_totsig - y_lagsig)/y_lagsig
    y_max = np.max(y_deltasig)
    gs = gridspec.GridSpec(2,1, height_ratios=[4,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.plot(evolve_x[::my_slice], y_lagsig[::my_slice], 'b.', label='lag')
    ax1.plot(evolve_x[::my_slice], y_totsig[::my_slice], 'r.', label='tot')
    ax1.grid(True)
    ax1.set(ylabel="Z")
    ax1.set_title("Z evolve history")
    ax1.legend()
    ax2.grid(True)
    ax2.set(ylabel=r'$\frac{tot-lag}{lag}$%')
    ax2.plot(evolve_x[::my_slice],y_deltasig[::my_slice], 'k.')
    ax2.set(xlabel="#cell")
    ax2.set_ylim(-0.1*y_max,1.1*y_max)
    plt.savefig("Z_evolve_history_PairNegWgt{}_Delta{}_GroupEvt{}.png".format(negative_weight,delta,group_events))
    plt.clf()
