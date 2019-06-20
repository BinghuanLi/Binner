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
import ROOT
from ROOT import kBlack, kBlue, kRed, kCyan, kViolet, kGreen, kOrange, kGray, kPink, kTRUE
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TString, TFile, TH1, THStack, gROOT, TStyle, TAttFill, TLegend, TGraphAsymmErrors

# set mimimum bkg for significance calculation
minimum_bkg = 0.000001
negative_weight = "2DHist" # Sort/Dist
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
NXQ = 10
NYQ = 10

# load data
inputPath = "../data/"
outputPath = "../output/2DHist_NXQ{}_NYQ{}/".format(NXQ,NYQ)
plotPath = "../plots/2DHist_NXQ{}_NYQ{}_TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}/".format(NXQ,NYQ,tf,negative_weight,delta,group_events)
if not os.path.exists(plotPath):
    os.makedirs(plotPath)
if not os.path.exists(outputPath):
    os.makedirs(outputPath)
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
my_data = my_data.astype({"bin_region":np.int64})
data = my_data.sample(frac=tf,random_state=200)
df_test = my_data.drop(data.index)
data = data.reset_index(drop=True)
df_test = df_test.reset_index(drop=True)


# load signal points to generate 2D diagram bin edges
sig_points_x = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttV"].values
sig_points_y = data.ix[(data.target.values==1)]["mvaOutput_2lss_ttbar"].values
sig_wgt = data.ix[(data.target.values==1)]["totalWeight"].values

print ("calculate weighted quantiles")
XQ = tool.get_wtd_quantiles(sig_points_x,sig_wgt,NXQ,500,-1,1)
YQ = tool.get_wtd_quantiles(sig_points_y,sig_wgt,NYQ,500,-1,1)
print(XQ)
print(YQ)
print ("finish calculate weighted quantiles")

# calculate signals
my_s, my_sumw2_s, my_entrie_s, df_bin = tool.fill_2D_hist(sig_points_x,sig_points_y,sig_wgt,XQ,YQ)

# calculate b1 (ttV) 
ttV_x = data.ix[(data.target.values==2)]["mvaOutput_2lss_ttV"].values
ttV_y = data.ix[(data.target.values==2)]["mvaOutput_2lss_ttbar"].values
ttV_wgt = data.ix[(data.target.values==2)]["totalWeight"].values
my_b1, my_sumw2_b1, my_entrie_b1, df_none = tool.fill_2D_hist(ttV_x,ttV_y,ttV_wgt,XQ,YQ)

# calculate b2 (ttJ) 
ttJ_x = data.ix[(data.target.values==3)]["mvaOutput_2lss_ttV"].values
ttJ_y = data.ix[(data.target.values==3)]["mvaOutput_2lss_ttbar"].values
ttJ_wgt = data.ix[(data.target.values==3)]["totalWeight"].values
my_b2, my_sumw2_b2, my_entrie_b2, df_none = tool.fill_2D_hist(ttJ_x,ttJ_y,ttJ_wgt,XQ,YQ)

print ("finish assigning bkgs")
'''
if doPlot:
    print (" plot initial 2D distributions ")
    
    plt.hist2d(sig_points_x,sig_points_y,bins=[XQ,YQ],weights=sig_wgt,cmap='rainbow')
    plt.colorbar()
    plt.title("Signal ttH")
    plt.xlabel("mvaOutput_2lss_ttV")
    plt.ylabel("mvaOutput_2lss_ttbar")
    plt.savefig("{}TrainFrac{}_2DHisto_ttH_initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()
    
    
    plt.hist2d(ttV_x,ttV_y,bins=[XQ,YQ],weights=ttV_wgt,cmap='rainbow')
    plt.colorbar()
    plt.title("Bkg1 ttV")
    plt.xlabel("mvaOutput_2lss_ttV")
    plt.ylabel("mvaOutput_2lss_ttbar")
    plt.savefig("{}TrainFrac{}_2DHisto_ttV_initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()

    plt.hist2d(ttJ_x,ttJ_y,bins=[XQ,YQ],weights=ttJ_wgt,cmap='rainbow')
    plt.colorbar()
    plt.title("Bkg2 ttJ")
    plt.xlabel("mvaOutput_2lss_ttV")
    plt.ylabel("mvaOutput_2lss_ttbar")
    plt.savefig("{}TrainFrac{}_2DHisto_ttJ_initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()

    print (" finish plot initial 2D histogram ")
'''

my_bin_region = np.arange(len(my_s)) 
my_err_sum_b = np.zeros(len(my_s))
my_significance = np.zeros(len(my_s))
#print(my_err_sum_b)
print ("calculate significance")
vfunc = np.vectorize(tool.get_significance)
my_err_sum_b, my_significance = vfunc(my_s, my_b1, my_b2, err_b1, err_b2, minimum_bkg)
print ("finish calculate significance")
df_hist =  pd.DataFrame({"bin_region":my_bin_region,"x_center":df_bin.x_center.values,"x_lowedge":df_bin.x_lowedge.values,"x_highedge":df_bin.x_highedge.values,"y_center":df_bin.y_center.values,"y_lowedge":df_bin.y_lowedge.values,"y_highedge":df_bin.y_highedge.values,"s":my_s,"b1":my_b1,"b2":my_b2,"sumw2_s":my_sumw2_s,"sumw2_b1":my_sumw2_b1,"sumw2_b2":my_sumw2_b2,"entrie_s":my_entrie_s,"entrie_b1":my_entrie_b1,"entrie_b2":my_entrie_b2,"error":my_err_sum_b,"significance":my_significance})


# sort by significance
df_hist.sort_values(by="significance", ascending=False, inplace=True)
df_hist = df_hist.reset_index(drop=True)

file_log.write("finish calculate siginificance, after sorting df_hist is {}\n".format(df_hist))

if doPlot:
    # plot the original significance
    print (" plot original significance ")
    print (df_hist["significance"].values)
    print (df_hist["significance"].values[0])
    values1, bins, _ = plt.hist(x=df_hist["significance"].values, bins = 10, alpha =0.5, range=(0., 1.01*df_hist["significance"].values[0]))
    plt.xlabel("Z")
    plt.ylabel("Entries")
    plt.title("Z distribution")
    plt.savefig("{}TrainFrac{}_Significance_Initial_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
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
n_cells = len(df_hist.index)
# set cost function
cost_fom = 0
sum_Zsquare =  (df_hist["significance"]**2).sum()
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
for index, row in df_hist.iterrows():
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

vor_cell = vor_cell+1
df_hist["total_sig"] = total_Z
df_hist["lag_sig"] = lag_Z
df_hist["update"] = update_Z
df_hist["bin_label"] = vor_cell


df_hist.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_df_hist.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))
print("finish merge algorithm")
file_log.write("finish merge algorithm df_hist is saved into {}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_df_hist.csv \n".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))

# save map
#print(data)
df_map =  df_hist[["x_center","x_lowedge","x_highedge","y_center","y_lowedge","y_highedge","bin_region","bin_label"]] 
print("df_map is created ")
df_map.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_map.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))
h2 = tool.map_to_2D_hist(df_map,XQ,YQ)
f1 = ROOT.TFile.Open("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_map.root".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr),"recreate")
f1.cd()
h2.Write()
f1.Close()
print("finish save map")

# ksTest
# create train df
df_hist_train = df_hist.groupby('bin_label').sum()[['s','sumw2_s','b1','sumw2_b1','b2','sumw2_b2','significance']].reset_index()
#df_hist_train = df_hist.groupby('bin_label').sum()[['s','sumw2_s','b1','sumw2_b1','b2','sumw2_b2']]
_, df_hist_train['significance'] = vfunc(df_hist_train['s'].values, df_hist_train['b1'].values, df_hist_train['b2'].values, err_b1, err_b2, minimum_bkg)
#_, my_train_significance = vfunc(df_hist_train['s'].values, df_hist_train['b1'].values, df_hist_train['b2'].values, err_b1, err_b2, minimum_bkg)
#df_hist_train['significance']=my_train_significance
#print(df_hist_train)
#print(df_hist_train['significance'])
print("finish recalculate significance")
file_log.write(" overtrain test: \n df_hist_train is - \n {} \n".format(df_hist_train))
print("finish load train data")
# create test df
# apply the map:
bin_data = tool.apply_bin_map_root(df_test,"{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_map.root".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr),"h2")
bin_data.to_csv("{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_test.csv".format(outputPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr))
print("finish load test data")
# to_table
test_values = ["bin_region","totalWeight","error","entries","sumw2"]
test_index = ["bin_label","target"]
#print(df)
test_data_table = tool.to_table(bin_data,test_values,test_index)
test_bin_label = test_data_table.index.get_level_values(0)[test_data_table.index.get_level_values(1)==1].values
test_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["totalWeight"].values
test_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["totalWeight"].values
test_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["totalWeight"].values
test_sumw2_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["sumw2"].values
test_sumw2_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["sumw2"].values
test_sumw2_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["sumw2"].values
test_entrie_s = test_data_table.loc[test_data_table.index.get_level_values(1)==1]["entries"].values
test_entrie_b1 = test_data_table.loc[test_data_table.index.get_level_values(1)==2]["entries"].values
test_entrie_b2 = test_data_table.loc[test_data_table.index.get_level_values(1)==3]["entries"].values
test_significance = np.zeros(len(test_bin_label))
_, test_significance = vfunc(test_s, test_b1, test_b2, err_b1, err_b2, minimum_bkg)
df_hist_test =  pd.DataFrame({"bin_label":test_bin_label,"s":test_s,"b1":test_b1,"b2":test_b2,"sumw2_s":test_sumw2_s,"sumw2_b1":test_sumw2_b1,"sumw2_b2":test_sumw2_b2,"significance":test_significance})
#df_hist_test.set_index('bin_label',inplace=True)
file_log.write(" overtrain test: \n df_hist_test is - \n {}\n".format(df_hist_test))

print(" finish organizing train and test data")
# make plots
#tool.make_compare_plot(df_hist_train, df_hist_test, ["s","b1","b2","significance"], "{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_TrainVsTest_{}".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
#tool.make_compare_plot(df_hist_train, df_hist_test, ["s"], "{}TrainFrac{}_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_TrainVsTest_{}".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))

print(" finish plot train/test comparison and prepare to plot final 2D bin region ")

'''
if doPlot :
    print (" plot final merged 2D diagrams ")
    
    
    c = TCanvas("c", "canvas", 800, 800)
    c.cd()
    h2.GetXaxis().SetTitle("mvaOutput_2lss_ttV")
    h2.GetYaxis().SetTitle("mvaOutput_2lss_ttbar")
    h2.SetStats(0)
    h2.Draw("COLZ") 
    c.SaveAs("{}TrainFrac{}_Final2DHist_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))

if doPlot: 
    # plot the evolving history
    print (" plot the evolving history ")
    my_slice = int(math.ceil(n_cells/100.))
    evolve_x = np.arange(n_cells)
    y_totsig = df_hist["total_sig"].values
    y_lagsig = df_hist["lag_sig"].values
    y_deltasig = 100.*(y_totsig - y_lagsig)/y_lagsig
    print("finish load data")
    y_max = np.max(y_deltasig[::my_slice])
    y_min = np.min(y_deltasig[::my_slice])
    gs = gridspec.GridSpec(2,1, height_ratios=[4,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    print("finish create axes")
    #ax1.plot(evolve_x[::my_slice], y_lagsig[::my_slice], 'b.', label='lag')
    #ax1.plot(evolve_x[::my_slice], y_totsig[::my_slice], 'r.', label='tot')
    #ax1.plot(evolve_x, y_lagsig, 'b.', label='lag')
    #ax1.plot(evolve_x, y_totsig, 'r.', label='tot')
    ax1.plot([1.,2.], [1.,1.], 'b.', label='lag')
    ax1.plot([1.,2.], [2.,2.], 'r.', label='tot')
    ax1.grid(True)
    ax1.set(ylabel="Z")
    ax1.set_title("Z evolve history")
    ax1.legend()
    print("finish ax1")
    ax2.grid(True)
    ax2.set(ylabel=r'$\frac{tot-lag}{lag}$%')
    #ax2.plot(evolve_x[::my_slice],y_deltasig[::my_slice], 'k.')
    ax2.plot([1.,2.],[1.2,1.4], 'k.')
    print("finish plot ax2")
    ax2.set(xlabel="#cell")
    #ax2.set_ylim(y_min,1.1*y_max)
    ax2.set_yscale('log')
    print("finish scale ax2")
    plt.savefig("{}TrainFrac{}_2DHist_Z_evolve_history_PairNegWgt{}_Delta{}_GroupEvt{}_useErr{}_minMC{}_maxErr{}_{}.png".format(plotPath,tf,negative_weight,delta,group_events,useError,minMC,maxErr,date))
    plt.close()
    print("finish close plt")
'''
file_log.close()
print("finish close file log")
