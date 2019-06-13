# coding: utf-8
from datetime import datetime
import os
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
delta = 1 # load event every delta instance, so the events loaded would be N/delta
group_events = -1 # group events so each initial cell contains group_events instance, the initial #cell = (N/delta)/group_events
doPlot = True
load_csv = True
err_b1 = 0.15 
err_b2 = 0.25 
tf = 0.5
useError = False
maxErr = 999
minMC = -1
date=datetime.today().strftime('%Y%m%d')

# load data
inputPath = "../data/"
outputPath = "../output/QT100/"
plotPath = "../plots/QT100_TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}/".format(tf,negative_weight,delta,group_events)
if not os.path.exists(plotPath):
    os.makedirs(plotPath)
variables = ["mvaOutput_2lss_ttV","mvaOutput_2lss_ttbar"]
# save log file
file_log = open(("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}.log".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr)),"w")

if not load_csv:
    my_data=load_data_2017(inputPath, variables, "passGenMatchCut==1", skip = delta ,nEvt=group_events, pair_negwgt = negative_weight ) 
    my_data.to_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}.csv".format(inputPath,negative_weight,delta,group_events,useError,minMC,maxErr))
else:
    try: my_data = pd.read_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}.csv".format(inputPath,negative_weight,delta,group_events,useError,minMC,maxErr))
    except:
        file_log.write( "{}input_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}.csv deosn't exist\n ".format(inputPath,negative_weight,delta,group_events,useError,minMC,maxErr))
        file_log.write( "run load_data_2017 \n")
        my_data=load_data_2017(inputPath, variables, "passGenMatchCut==1", skip = delta ,nEvt=group_events, pair_negwgt = negative_weight ) 
        my_data.to_csv("{}input_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}.csv".format(inputPath,negative_weight,delta,group_events,useError,minMC,maxErr))

print("finish load lata my_data")
file_log.write(str(my_data)+"\n")


# train test split
my_data = my_data.astype({"vor_point":np.int64,"vor_region":np.int64})
data = my_data.sample(frac=tf,random_state=200)
df_test = my_data.drop(data.index)
data = data.reset_index(drop=True)
df_test = df_test.reset_index(drop=True)


# load signal points to generate voronoi diagram
sig_points_x = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttV"]
sig_points_y = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttbar"]

Points=np.vstack((sig_points_x,sig_points_y)).T

Vor = Voronoi(Points)
#print("print voronoi vertices")
#print(Vor.vertices)
#print("print voronoi regions")
#print(Vor.regions)

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
#print("region idx is {}".format(region_idx))
mask = data.target ==2
data.loc[mask,"vor_region"]= region_idx
data.loc[mask,"vor_point"]= point_idx

point_idx, region_idx = tool.find_vor_region(ttJ_points,Vor)
#print("region idx is {}".format(region_idx))
mask = data.target ==3
data.loc[mask,"vor_region"]= region_idx
data.loc[mask,"vor_point"]= point_idx

# sort data by vor_region
print ("finish assigning bkgs")
file_log.write("finish assining bkgs, data is {} \n".format(data))

if doPlot and len(Vor.regions) < 1000 :
#if doPlot:
    # plot initial Voronoi diagrams only if vor.regions < 1000
    print (" plot initial Voronoi diagrams ")
    plt = tool.plot_vor_2d(Points, Vor.vertices, Vor.regions, Vor.ridge_vertices, Vor.ridge_points) 
    plt.xlabel("mvaOutput_2lss_ttV")
    plt.ylabel("mvaOutput_2lss_ttbar")
    plt.plot(ttV_x,ttV_y,'r.', label='ttV')
    plt.plot(ttJ_x,ttJ_y,'k.', label='ttJ')
    plt.legend()
    plt.title("initial voronoi diagrams")
    #plt.show()
    plt.savefig("{}TrainFrac{}_Voronoi_initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()

    print (" finish plot initial Voronoi diagrams ")

#column_list = ["vor_point","vor_region","s","b1","b2","b1_err","b2_err","error","significance"]
my_values = ["vor_point","totalWeight","error","entries","sumw2"]
my_index = ["vor_region","target"]
#print(df)
data_table = tool.to_table(data,my_values,my_index)
my_vor_region = data_table.index.get_level_values(0)[data_table.index.get_level_values(1)==1].values
#print (my_vor_region)
my_vor_point = data_table.loc[data_table.index.get_level_values(1)==1]["vor_point"].values
#print (my_vor_point)
my_s = data_table.loc[data_table.index.get_level_values(1)==1]["totalWeight"].values
my_b1 = data_table.loc[data_table.index.get_level_values(1)==2]["totalWeight"].values
my_b2 = data_table.loc[data_table.index.get_level_values(1)==3]["totalWeight"].values
my_sumw2_s = data_table.loc[data_table.index.get_level_values(1)==1]["sumw2"].values
my_sumw2_b1 = data_table.loc[data_table.index.get_level_values(1)==2]["sumw2"].values
my_sumw2_b2 = data_table.loc[data_table.index.get_level_values(1)==3]["sumw2"].values
my_entrie_s = data_table.loc[data_table.index.get_level_values(1)==1]["entries"].values
my_entrie_b1 = data_table.loc[data_table.index.get_level_values(1)==2]["entries"].values
my_entrie_b2 = data_table.loc[data_table.index.get_level_values(1)==3]["entries"].values
my_err_sum_b = np.zeros(len(my_vor_region))
my_significance = np.zeros(len(my_vor_region))
#print(my_err_sum_b)
print ("calculate significance")
vfunc = np.vectorize(tool.get_significance)
my_err_sum_b, my_significance = vfunc(my_s, my_b1, my_b2, err_b1, err_b2, minimum_bkg)

print ("finish calculate significance")

df_vor =  pd.DataFrame({"vor_point":my_vor_point,"vor_region":my_vor_region,"s":my_s,"b1":my_b1,"b2":my_b2,"sumw2_s":my_sumw2_s,"sumw2_b1":my_sumw2_b1,"sumw2_b2":my_sumw2_b2,"entrie_s":my_entrie_s,"entrie_b1":my_entrie_b1,"entrie_b2":my_entrie_b2,"error":my_err_sum_b,"significance":my_significance})


# sort by significance
df_vor.sort_values(by="significance", ascending=False, inplace=True)
df_vor = df_vor.reset_index(drop=True)

file_log.write("finish calculate siginificance, after sorting df_vor is {}\n".format(df_vor))

if doPlot:
    # plot the original significance
    print (" plot original significance ")
    values1, bins, _ = plt.hist(df_vor["significance"].values, bins = 10, alpha =0.5, range=(0., 1.01*df_vor["significance"].values[0]) , density = True)
    plt.xlabel("Z")
    plt.ylabel("%")
    plt.title("Z distribution")
    plt.savefig("{}TrainFrac{}_Significance_Initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()

# get and add qunatiles to df_vor
vor_qt = pd.qcut(range(len(df_vor)),100,labels=False)
df_vor["vor_qt"] = vor_qt
# grouping by quantiles
df_quantile = df_vor[['s','b1','b2','sumw2_s','sumw2_b1','sumw2_b2','entrie_s','entrie_b1','entrie_b2','vor_qt']]
print(df_quantile)
df_quantile = df_quantile.groupby('vor_qt').sum().reset_index()
print(df_quantile)

print ("calculate qt significance")
my_err_sum_b, qt_significance = vfunc(df_quantile.s.values, df_quantile.b1.values, df_quantile.b2.values, err_b1, err_b2, minimum_bkg)
df_quantile["significance"] = qt_significance
df_quantile.sort_values(by="significance", ascending=False, inplace=True)
df_quantile = df_quantile.reset_index(drop=True)

print ("finish calculate qt significance")

if doPlot:
    # plot the qt significance
    print (" plot qt significance ")
    values1, bins, _ = plt.hist(df_quantile["significance"].values, bins = 10, alpha =0.5, range=(0., 1.01*df_quantile["significance"].values[0]) , density = True)
    plt.xlabel("Z")
    plt.ylabel("%")
    plt.title("Z distribution")
    plt.savefig("{}TrainFrac{}_QT_Significance_Initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()

# bin merge algorithm
s_trail =0
b1_trail = 0
b2_trail = 0
sumw2_b1_trail = 0
sumw2_b2_trail = 0
entrie_b1_trail = 0
entrie_b2_trail = 0
Z_trail = 0
n_cells = len(df_quantile.index)
# set cost function
cost_fom = 0
sum_Zsquare =  (df_quantile["significance"]**2).sum()
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
for index, row in df_quantile.iterrows():
    s = row["s"]
    b1 = row["b1"]
    b2 = row["b2"]
    sumw2_b1 = row["sumw2_b1"]
    sumw2_b2 = row["sumw2_b2"]
    entrie_b1 = row["entrie_b1"]
    entrie_b2 = row["entrie_b2"]
    Z = row["significance"]
    lag_Z_square = sum_Zsquare
    lag_delta_Z = sum_deltaZ
    err_sum_b , significance = tool.get_significance((s_trail+s), (b1_trail+b1), (b2_trail+b2), err_b1, err_b2, minimum_bkg)
    #print("s_trail:{}, b1_trail:{}, b2_trail:{}, s:{}, b1:{}, b2:{}, Cost:{}, Significance: {}, Sum_Zsqure: {}".format(s_trail, b1_trail, b2_trail, s, b1, b2, cost,(significance**2),(Z**2+Z_trail**2)))
    B = b1+b1_trail+b2+b2_trail
    sumw2_B = sumw2_b1+sumw2_b1_trail+sumw2_b2+sumw2_b2_trail
    entrie_B = entrie_b1+entrie_b1_trail+entrie_b2+entrie_b2_trail
    enough_bkg = tool.have_enough_stat(B, sumw2_B, entrie_B, minMC,maxErr, useErr = useError)
    cost = tool.cost_fun(significance, Z_trail, Z, cost_fom)
    if cost >= 0. or not enough_bkg:
        # merge the bins
        sum_Zsquare += cost
        sum_deltaZ += cost
        s_trail += s
        b1_trail += b1
        b2_trail += b2
        sumw2_b1_trail += sumw2_b1
        sumw2_b2_trail += sumw2_b2
        entrie_b1_trail += entrie_b1
        entrie_b2_trail += entrie_b2
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
        sumw2_b1_trail = sumw2_b1
        sumw2_b2_trail = sumw2_b2
        entrie_b1_trail = entrie_b1
        entrie_b2_trail = entrie_b2
        Z_trail = Z
        if cost_fom ==0:
            lag_Z[index] = math.sqrt(lag_Z_square)
            total_Z[index] = math.sqrt(sum_Zsquare)
        else:
            lag_Z[index] = sum_deltaZ
            total_Z[index] = sum_deltaZ
        vor_count +=1
    vor_cell[index] = vor_count 

df_quantile["total_sig"] = total_Z
df_quantile["lag_sig"] = lag_Z
df_quantile["update"] = update_Z
df_quantile["vor_label"] = vor_cell

print(df_quantile)

# now save it to the DataFrame
all_qt_labels = dict(zip(df_quantile["vor_qt"].values,vor_cell))
print( " vor_qt : vor_cell dictionary ")
print(df_quantile["vor_qt"].values)
print(vor_cell)
print(all_qt_labels)
df_vor['vor_label'] = df_vor['vor_qt'].map(all_qt_labels)

df_vor.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_df_vor.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))
print("finish merge algorithm")
file_log.write("finish merge algorithm df_vor is saved into {}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_df_vor.csv \n".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))

# dictionary to save {label:[p1,p2,..., pN]}
label_points = df_vor.groupby('vor_label')['vor_point'].apply(list).to_dict() 

regions, vertices = voronoi_finite_polygons_2d(Vor.points, Vor.vertices, Vor.regions, Vor.ridge_vertices, Vor.ridge_points, Vor.point_region)

# save map
#print(data)
df_map, drop_x, drop_y = tool.save_vor_map(data.ix[(data.target.values==1)], regions, vertices, label_points)
df_map = df_map.reset_index(drop=True)
df_map.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_map.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr), columns=["mvaOutput_2lss_ttV","mvaOutput_2lss_ttbar","vor_label"])
    
# ksTest
# create train df
df_vor_train = df_vor.groupby('vor_label').sum()[['s','sumw2_s','b1','sumw2_b1','b2','sumw2_b2','significance']]
_, df_vor_train['significance'] = vfunc(df_vor_train['s'], df_vor_train['b1'], df_vor_train['b2'], err_b1, err_b2, minimum_bkg)
file_log.write(" overtrain test: \n df_vor_train is - \n {} \n".format(df_vor_train))
print("finish load train data")
# create test df
# apply the map:
vor_data = tool.apply_vor_map(df_test,df_map)
vor_data.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_test.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))
print("finish load test data")
# to_table
test_values = ["vor_point","totalWeight","error","entries","sumw2"]
test_index = ["vor_label","target"]
#print(df)
test_data_table = tool.to_table(vor_data,test_values,test_index)
test_vor_label = test_data_table.index.get_level_values(0)[test_data_table.index.get_level_values(1)==1].values
test_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["totalWeight"].values
test_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["totalWeight"].values
test_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["totalWeight"].values
test_sumw2_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["sumw2"].values
test_sumw2_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["sumw2"].values
test_sumw2_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["sumw2"].values
test_entrie_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["entries"].values
test_entrie_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["entries"].values
test_entrie_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["entries"].values
test_significance = np.zeros(len(test_vor_label))
_, test_significance = vfunc(test_s, test_b1, test_b2, err_b1, err_b2, minimum_bkg)
df_vor_test =  pd.DataFrame({"vor_label":test_vor_label,"s":test_s,"b1":test_b1,"b2":test_b2,"sumw2_s":test_sumw2_s,"sumw2_b1":test_sumw2_b1,"sumw2_b2":test_sumw2_b2,"significance":test_significance})
df_vor_test.set_index('vor_label',inplace=True)
file_log.write(" overtrain test: \n df_vor_test is - \n {}\n".format(df_vor_test))

print(" finish organizing train and test data")
# make plots
tool.make_compare_plot(df_vor_train, df_vor_test, ["s","b1","b2","significance"], "{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_TrainVsTest_{}".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))

print(" finish plot train/test comparison and prepare to plot final voronoi region ")

if doPlot and len(Vor.regions) < 5000 :
#if doPlot :
    # plot merged Voronoi diagrams with color
    print (" plot final merged Voronoi diagrams ")
    

    #print ("regions")
    #print (regions)
    #print ("vertices")
    #print (vertices)

    final_label_size = df_vor["vor_label"].nunique()
    final_labels=df_vor.vor_label.unique()
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
    #plt.plot(drop_x,drop_y,'k.', label='drop')
    plt.legend()

    #plt.show()
    plt.savefig("{}TrainFrac{}_ColoredVoronoi_Final_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()
    print(" finish plot final voronoi region and prepare to plot evolving history ")

if doPlot: 
    # plot the evolving history
    print (" plot the evolving history ")
    my_slice = int(math.ceil(n_cells/100.))
    evolve_x = np.arange(n_cells)
    y_totsig = df_quantile["total_sig"].values
    y_lagsig = df_quantile["lag_sig"].values
    y_deltasig = 100.*(y_totsig - y_lagsig)/y_lagsig
    y_max = np.max(y_deltasig[::my_slice])
    y_min = np.min(y_deltasig[::my_slice])
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
    ax2.set_ylim(min(-0.1*y_max, 1.1*y_min),1.1*y_max)
    plt.savefig("{}TrainFrac{}_QT_Z_evolve_history_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.clf()


file_log.close()
