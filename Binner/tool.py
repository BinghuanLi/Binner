""" Tools """

import ROOT
import numpy as np
import pandas as pd
from root_numpy import root2array, tree2array 
from scipy.spatial import KDTree, cKDTree
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import math
from ROOT import RooStats
from scipy import stats
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib import gridspec
from shapely import geometry
from itertools import compress
import shapely
import geopandas
from collections import deque

def bin_events(df, variables, nEvt =-1):
    '''
        sort the events by distance and then group every nEvt events
        input : dataFrame, variables, nCell
        output : grouped dataFrame 
    '''
    if nEvt <=1:
        return df
    # sort by sum of two BDT values
    df = df.ix[np.argsort(df.mvaOutput_2lss_ttV + df.mvaOutput_2lss_ttbar).values]
    df.reset_index(drop=True, inplace = True)
    #print(df)
    my_list_column = variables + ['totalWeight','entries']
    my_list_agg = ['mean']*len(variables)
    my_dict = dict(zip(my_list_column,my_list_agg))
    my_dict['totalWeight']='sum'
    my_dict['entries']='sum'
    #print(my_dict)
    df = df.groupby(df.index / nEvt).agg(my_dict)
    return df

def pair_negative_distance(df):
    '''
        pair non positive weighted events to closest point 
        super slow
    '''
    step = math.sqrt(2./(len(df.index)+1))
    print(df)
    print(step)
    # get the indices of non positive weights
    neg_indices = df.index[df['entries']<=0].tolist()
    print(neg_indices)
    # the indices of rows to be dropped
    drop_indices = []
    # get points 
    Points = np.vstack((df["mvaOutput_2lss_ttV"],df["mvaOutput_2lss_ttbar"])).T
    for i in neg_indices:
        print (i)
        #print( 'dropped indices')
        #print( drop_indices)
        if i in drop_indices: continue
        tree = KDTree(Points)
        coordinate = [df.at[i,'mvaOutput_2lss_ttV'],df.at[i,'mvaOutput_2lss_ttbar']]
        dist=0
        while dist < 10*step and (df.at[i,'entries'] <=0 or df.at[i,'totalWeight'] <=0):
            dist += step
            point_idx = tree.query_ball_point(coordinate,dist)
            #print (point_idx)
            if df.loc[point_idx,'entries'].sum() >0 and df.loc[point_idx,'totalWeight'].sum() >0:
                df.at[i,'totalWeight'] = df.loc[point_idx,'totalWeight'].sum()
                df.at[i,'entries'] = df.loc[point_idx,'entries'].sum()
                df.at[i,'mvaOutput_2lss_ttV'] = df.loc[point_idx,'mvaOutput_2lss_ttV'].mean()
                df.at[i,'mvaOutput_2lss_ttbar'] = df.loc[point_idx,'mvaOutput_2lss_ttbar'].mean()
                point_idx.remove(i)
                drop_indices += point_idx
    
    df.drop(drop_indices, inplace = True)
    df.reset_index(drop=True, inplace = True)
    neg_indices = df.index[df['entries']<=0].tolist()
    if (len(neg_indices)>0):
        print("NOTICE: still have negative weighted points, please change the step and check the events ")
    print(df)
    return df

def pair_negative(df):
    '''
        pair non positive weighted events to nearby points 
    '''
    # sort by sum of two BDT values
    df = df.ix[np.argsort(df.mvaOutput_2lss_ttV + df.mvaOutput_2lss_ttbar).values]
    df.reset_index(drop=True, inplace = True)
    print (df)
    # the indices of rows to be dropped
    drop_indices = []
    # get the indices of non positive weights
    neg_indices = df.index[df['entries']<=0].tolist()
    # get the looped indices of non positive weights
    loop_indices = [] 
    # loop forward
    #print(neg_indices)
    for i in neg_indices:
        if i in drop_indices: continue
        count = 1
        while (count+i) <=  df.index.max() :
            df.at[i,'entries'] = df.at[i,'entries']+df.at[count+i,'entries']
            df.at[i,'totalWeight'] = df.at[i,'totalWeight']+ df.at[count+i,'totalWeight']
            df.at[i,'mvaOutput_2lss_ttV'] = df.at[i,'mvaOutput_2lss_ttV']+ df.at[count+i,'mvaOutput_2lss_ttV']
            df.at[i,'mvaOutput_2lss_ttbar'] = df.at[i,'mvaOutput_2lss_ttbar']+ df.at[count+i,'mvaOutput_2lss_ttbar']
            drop_indices.append(count+i)
            count +=1
            if df.at[i,'entries'] > 0 and df.at[i,'totalWeight'] > 0:
                df.at[i,'mvaOutput_2lss_ttV'] = df.at[i,'mvaOutput_2lss_ttV']/count 
                df.at[i,'mvaOutput_2lss_ttbar'] = df.at[i,'mvaOutput_2lss_ttbar']/count 
                loop_indices.append(i)
                break
    #print(drop_indices)
    #print(neg_indices)
    #print(loop_indices)
    df.drop(drop_indices, inplace = True)
    df.reset_index(drop=True, inplace = True)
    
    # loop backward
    if len(neg_indices)> len(loop_indices):
        print ( 'loop backward')
        # the indices of rows to be dropped
        drop_indices = []
        # get the indices of non positive weights
        neg_indices = df.index[df['entries']<=0].tolist()
        # get the looped indices of non positive weights
        loop_indices = [] 
        print(neg_indices)
        for i in neg_indices:
            if i in drop_indices: continue
            count = 1
            while (i-count) >=  df.index.min():
                df.at[i,'entries'] = df.at[i,'entries']+df.at[i-count,'entries']
                df.at[i,'totalWeight'] = df.at[i,'totalWeight']+ df.at[i-count,'totalWeight']
                df.at[i,'mvaOutput_2lss_ttV'] = df.at[i,'mvaOutput_2lss_ttV']+ df.at[i-count,'mvaOutput_2lss_ttV']
                df.at[i,'mvaOutput_2lss_ttbar'] = df.at[i,'mvaOutput_2lss_ttbar']+ df.at[i-count,'mvaOutput_2lss_ttbar']
                drop_indices.append(i-count)
                count +=1
                if df.at[i,'entries'] > 0 and df.at[i,'totalWeight'] > 0: 
                    df.at[i,'mvaOutput_2lss_ttV'] = df.at[i,'mvaOutput_2lss_ttV']/count 
                    df.at[i,'mvaOutput_2lss_ttbar'] = df.at[i,'mvaOutput_2lss_ttbar']/count 
                    loop_indices.append(i)
                    break
        print(drop_indices)
        print(neg_indices)
        df.drop(drop_indices, inplace = True)
        df.reset_index(drop=True, inplace = True)
    
    # drop the remaining neg_indices
    if len(neg_indices)> len(loop_indices):
        print( 'I have loop over forward and backward but neg_indices still exists: ')
        # get the indices of non positive weights
        neg_indices = df.index[df['entries']<=0].tolist()
        df.drop(neg_indices, inplace = True)
        df.reset_index(drop=True, inplace = True)
        print(drop_indices)
        print(neg_indices)
    
    
    print (df)
    return df


#load_data_2017(inputPath, variables, False)  # select all jets
def load_data_2017(inputPath,variables,criteria, skip=1, nEvt=-1, pair_negwgt = "Sort"):
    print variables
    my_cols_list=variables+['proces', 'key', 'target','totalWeight','entries','error','vor_point','vor_region']
    data = pd.DataFrame(columns=my_cols_list)
    keys=['ttHnobb_SigRegion','ttWJets_SigRegion',"ttZJets_SigRegion","ttJets_SigRegion"]
    print keys
    for key in keys :
        print key
        if 'ttH' in key or 'TTH' in key:
                sampleName='ttH'
                target=1
                error=0.
        if 'ttW' in key or 'TTW' in key or 'ttZ' in key or 'TTZ' in key:
                sampleName='ttV'
                target=2
                error=0.15
        if 'ttJets' in key or 'TTJets' in key:
                sampleName='ttJ'
                target=3
                error=0.25
        
        inputTree = 'mvaTree'
        try: tfile = ROOT.TFile(inputPath+"/"+key+".root")
        except : 
            print " file "+ inputPath+"/"+key+".root deosn't exits "
            continue
        try: tree = tfile.Get(inputTree)
        except : 
            print inputTree + " deosn't exists in " + inputPath+"/"+key+".root"
            continue
        if tree is not None :
            try:
                chunk_arr = tree2array(tree=tree, selection=criteria, step=skip)
            except : continue
            else :
                #print (chunk_arr)
                chunk_df = pd.DataFrame(chunk_arr, columns=variables)
                chunk_df['totalWeight']=chunk_arr['EventWeight']
                chunk_df['entries']=np.sign(chunk_arr['EventWeight'])
                if pair_negwgt == "Sort":
                    chunk_df=pair_negative(chunk_df) 
                elif pair_negwgt == "Dist":
                    chunk_df=pair_negative_distance(chunk_df) 
                else: 
                    print ("WARNING negative weight is not paired ")
                chunk_df= bin_events(chunk_df, variables, nEvt)
                #print (chunk_df)
                #print (chunk_df.columns.tolist())
                #print( "sampleName "+ sampleName)
                #print( "key "+ key)
                #print( "target "+ str(target))
                chunk_df['proces']=sampleName
                chunk_df['key']=key
                chunk_df['target']=target
                chunk_df['vor_region']=-1
                chunk_df['vor_point']=-1
                chunk_df['error']=error
                print(chunk_df)
                print(data)
                data=data.append(chunk_df, ignore_index=True, sort = True)
        tfile.Close()
        if len(data) == 0 : continue
        nS = len(data.ix[(data.target.values == 1) & (data.key.values==key) ])
        nB = len(data.ix[(data.target.values > 1) & (data.key.values==key) ])
        print key,"length of sig, bkg: ", nS, nB , data.ix[ (data.key.values==key)]["totalWeight"].sum()
        nNW = len(data.ix[(data["totalWeight"].values < 0) & (data.key.values==key) ])
        print key, "events with -ve weights", nNW
    print (data.columns.values.tolist())
    n = len(data)
    nS = len(data.ix[data.target.values == 1])
    nB1 = len(data.ix[data.target.values ==2 ])
    nB2 = len(data.ix[data.target.values ==3 ])
    print " length of data, sig, bkg1, bkg2: ", n, nS, nB1, nB2
    n = data["totalWeight"].sum()
    nS = data.ix[ (data.target.values==1)]["totalWeight"].sum() 
    nB1 = data.ix[ (data.target.values==2)]["totalWeight"].sum() 
    nB2 = data.ix[ (data.target.values==3)]["totalWeight"].sum() 
    print " yield of data, sig, bkg1, bkg2: ", n, nS, nB1, nB2
    return data


def find_vor_region(bkg_points, vor):
    """ find the voronoi region bkg_points belongs to """
    # get voronoi intinial points 
    Points = vor.points
    # get voronoi point_region
    Point_Region = vor.point_region
    # find the [x,y]'s closest neighbor 
    tree = KDTree(Points)
    # closest neighbor point index
    point_idx = tree.query(bkg_points)[1]
    # index vor region that closest neighbor point belongs to
    region_idx = Point_Region[point_idx]
    
    # print the neighbor point
    #print("neighbor point of {} is Point {}".format(bkg_points,Points[point_idx]))

    return point_idx, region_idx

def plot_vor_2d(points, vertices, regions, ridge_vertices, ridge_points ):
    ''' 
        the function to plot voronoi 
        reference : 
            https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#qhulltutorial
            https://gist.github.com/pv/8036995
        input:
            points:
                2D array
                The coordinates of points generate the voronoi diagram
            vertices: 
                2D array 
                The Voronoi vertices denote the set of points forming the polygonal edges of the Voronoi regions
            regions: 
                2D array
                These numbers indicate indices of the voronoi vertices making up the region, -1 indicates infinity
            ridge_vertices:
                2D array
                These numbers indicate indices of the Voronoi vertices making up the line segments. -1 is again a point at infinity
            ridge_points:
                2D array
                These numbers indicate indices of input points the voronoi ridges are perpendicular to 
        output:
            figure and plot
    '''
    # first plot the ppoints and the voronoi vertices
    plt.plot(points[:,0], points[:,1],'.',label='ttH')
    #plt.plot(vertices[:, 0], vertices[:, 1], '*')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)

    print("finish plot voronoi seed cells")
    # plot the finite line segments
    for simplex in ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex>=0):
            plt.plot(vertices[simplex, 0], vertices[simplex, 1], '-', color='orange', linewidth=2, alpha =0.6)

    print("finish plot the finite line segments")
    # plot the infinity
    center = points.mean(axis=0)
    for pointidx, simplex in zip(ridge_points, ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0] # finite end Voronoi vertex
            t = points[pointidx[1]] - points[pointidx[0]]  # tangent
            t = t / np.linalg.norm(t)
            n = np.array([-t[1], t[0]]) # normal
            midpoint = points[pointidx].mean(axis=0)
            far_point = vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
            plt.plot([vertices[i,0], far_point[0]],
                     [vertices[i,1], far_point[1]], '--', color='orange', linewidth=2, alpha =0.6)
    print("finish plot the infinite line segments")
    
    return plt

def eff_error(a, b, err_a, err_b, threshold, correlation=0):
    '''
        effective error 
    '''
    if (b+a) < threshold:
        err = 0
    else:
        err = math.sqrt((a*a*err_a*err_a + b*b*err_b*err_b + 2*correlation*a*b*err_a*err_b))/(b+a)
    return err

def get_significance(s, b1, b2, err_b1, err_b2, threshold):
    '''
        asimov significance
        https://root.cern.ch/doc/v614/RooStatsUtils_8cxx_source.html
        https://www.pp.rhul.ac.uk/~cowan/stat/medsig/medsigNote.pdf
        eq[10], eq[20]
    '''
    err_b1_b2 = eff_error(b1, b2, err_b1, err_b2, threshold)
    nB = max(b1+b2, threshold)
    significance = 0.
    if s >= 0:
        # please think and figure out a better way to treat negative weight 
        significance = RooStats.AsimovSignificance(s, nB, err_b1_b2)
    return err_b1_b2, significance

def cost_fun(z0, z1, z2, fom=0):
    '''
        calculate cost
        fom 0 : cost =  (z0**2) - (z1*2+z2**2)
        fom 1 : cost =  stats.norm.ppf(p_lead) -  stats.norm.ppf(p_trail)
    '''
    deltaZ = 0
    if fom==0:
        deltaZ = (z0**2) - (z1**2+z2**2)
    else:    
        p_trail = 1-(1-stats.norm.cdf(z1))*(1-stats.norm.cdf(z2))
        p_lead = stats.norm.cdf(z0)
        if (p_trail > 0.99999 and p_lead > 0.99999) or z1==0 : deltaZ = 0
        else: deltaZ = stats.norm.ppf(p_lead) -  stats.norm.ppf(p_trail) 
        print("z0,{},z1,{},z2,{},p_lead,{},p_trail,{},dZ,{}".format(z0,z1,z2,p_lead,p_trail,deltaZ))
     
    return deltaZ
    
def plot_polygon_collection(ax, geoms, values=None, colormap='gist_rainbow',  facecolor=None, edgecolor=None, alpha=0.5, linewidth=1.0, **kwargs):
    """ Plot a collection of Polygon geometries """
    patches = []
    for poly in geoms:
        a = np.asarray(poly.exterior)
        if poly.has_z:
            poly = shapely.geometry.Polygon(zip(*poly.exterior.xy))
        patches.append(Polygon(a))
    patches = PatchCollection(patches, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, alpha=alpha, **kwargs)

    if values is not None:
        patches.set_array(values)
        patches.set_cmap(colormap)
    ax.add_collection(patches, autolim=True)
    ax.autoscale_view()
    return patches

def count_ridges(x,y,z):
    return y.count(x) >=2 or len(set(x) & set(z))>=2

def my_list_count_ridges(x,y,z):            
    bool_list = [(y.count(item)>=2 or len(set(item) & set(z)) >=2)  for item in x]
    return all(bool_list)
 
def save_vor_map(data, new_regions, new_vertices, label_points):
    ''' 
        save the vor_map as dataFrame
        the idea is to:
            1. for each vor_label, create a label_ridges collections of {label:[ridges]}
            2. for each point in a given vor_label, if all the ridges within the plane of such point appears at least twice in the label_ridges collections, remove the point
        input : 
            data:
                DataFrame, contains the x,y coordinates of points
            new_regions:
                list of list of indicices of vor_region vertices, index of list is the vor_point number
            new_vertices:
                2D array contains the coordinate of vertices
            label_points:
                dictionary of {label:[p1,p2,..., pN]}

        output :
            df_map:
                DataFrame, contains minimum set of x,y,vor_label describing the voronoi map
    '''
    # create list save the index of vertices that is outside the 2D plane
    infinite_vertices = np.where(abs(new_vertices)>1)[0]
    
    # create the ridge pairs
    all_ridges = []
    for region in new_regions:
        old_region = region
        region = deque(region) 
        region.rotate(1)
        region=list(region)
        max_region = [ max(i,j) for i,j in zip(old_region,region) ]
        min_region = [ min(i,j) for i,j in zip(old_region,region) ]
        all_ridges.append([ list(a) for a in zip(max_region,min_region)])
   
    # loop over labels
    all_labels = []
    all_points = []
    for label, points in label_points.iteritems():  
        # create the label_ridge collection
        ridges = [all_ridges[x] for x in points]
        flat_ridges = [ item for sublist in ridges for item in sublist]
        # create the list of bool flag, true means the points needs to be save
        save_points= [ not my_list_count_ridges(ridge,flat_ridges,infinite_vertices) for ridge in ridges]
        # the list of points needs to be saved for this label
        my_points = list(compress(points,save_points))
        all_points = all_points + my_points
        all_labels = all_labels + [label] * len(my_points)

    # now save it to the DataFrame
    all_points_labels = dict(zip(all_points,all_labels))
    data['vor_label'] = data['vor_point'].map(all_points_labels)
    
    df_map = data.dropna(axis=0)
  
    df_na = data[data['vor_label'].isnull()]
    drop_x = df_na["mvaOutput_2lss_ttV"] 
    drop_y = df_na["mvaOutput_2lss_ttbar"] 
     
    return df_map, drop_x, drop_y
    

def to_table(data, my_values, my_index):
    '''
        convert dataFrame, so my_values is the function of my_list_agg and my_index is indices
    '''
    my_list_column = my_values 
    my_list_agg = ['mean']*len(my_list_column)
    my_dict = dict(zip(my_list_column,my_list_agg))
    my_dict['totalWeight']='sum'
    data_table = data.pivot_table(values=my_values, index=my_index, aggfunc = my_dict )
    data_table = data_table.unstack('target', fill_value=0).stack(level=1)
    return data_table


def apply_vor_map(data, vor_map):
    ''' 
        apply the vor_map to data
        input : 
            data:
                DataFrame, contains x,y and event weight
            vor_map:
                DataFrame, contains minimum set of x,y,vor_label describing the voronoi map

        output :
            vor_data:
                DataFrame, data with additional column vor_label
    '''
    # get voronoi points 
    Points = np.vstack((vor_map["mvaOutput_2lss_ttV"],vor_map["mvaOutput_2lss_ttbar"])).T
    # create the kDTree 
    tree = KDTree(Points)
    # create application points
    app_points = np.vstack((data["mvaOutput_2lss_ttV"],data["mvaOutput_2lss_ttbar"])).T
    # closest neighbor point index
    point_idx = tree.query(app_points)[1]
    print ( point_idx )
    # vor label that closest neighbor point belongs to
    labels = vor_map.loc[point_idx]["vor_label"].values
    vor_data = data
    vor_data["vor_label"] = labels

    return vor_data


def make_compare_plot(df_train, df_test, variables, figname, Norm=True, fig_sz=(10, 8)):
    '''
        OUTPUT: 
            No statistic test here because:
            - ks-Test is suitable only for continuous distribution 
            - Chi2 test is invalid when the obs or exp frequencies in each category are too small ( threshold 5)
        INPUTS:
        df_train - DataFrame voronoi train
        df_test - DataFrame voronoi application
        bins - number of bins for viz. Default 30.
        fig_sz - change to True in order to get larger outputs. Default False.
    
    '''
    
    x_train = df_train.index.values
    x_test = df_test.index.values
    for var in variables:
        y_train = df_train[var]
        y_test = df_test[var]
        if Norm:
            y_train = y_train/np.linalg.norm(y_train)
            y_test = y_test/np.linalg.norm(y_test)
    
        y_ratio=[(lambda x,y: x/y if y!=0 else(1 if x==0 else 0 ) )(a,b) for a,b in zip(y_train,y_test) ]


        gs = gridspec.GridSpec(2,1, height_ratios=[4,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax1.plot(x_train, y_train, 'b+', label='train')
        ax1.plot(x_test, y_test, 'rx', label='test')
        ax1.grid(True)
        ax1.set(ylabel="Normalized Unit")
        ax1.set_title(var)
        ax1.legend()
        ax2.grid(True)
        ax2.set(ylabel=r'$\frac{train}{test}$')
        ax2.plot(x_train,y_ratio, 'k.')
        ax2.set(xlabel="vor_label")
        ax2.set_ylim(0,2)
        plt.savefig("var-{}-{}".format(var,figname))
        plt.clf()
