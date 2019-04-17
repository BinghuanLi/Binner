""" Tools """

import ROOT
import numpy as np
import pandas as pd
from root_numpy import root2array, tree2array 
from scipy.spatial import KDTree
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import math
from ROOT import RooStats
from scipy import stats

#load_data_2017(inputPath, variables, False)  # select all jets
def load_data_2017(inputPath,variables,criteria):
    print variables
    my_cols_list=variables+['proces', 'key', 'target','totalWeight','error','vor_point','vor_region']
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
            try: chunk_arr = tree2array(tree=tree, selection=criteria, start=0, stop = 20)
            except : continue
            else :
                chunk_df = pd.DataFrame(chunk_arr, columns=variables)
                print (chunk_df.columns.tolist())
                #print( "sampleName "+ sampleName)
                #print( "key "+ key)
                #print( "target "+ str(target))
                chunk_df['proces']=sampleName
                chunk_df['key']=key
                chunk_df['target']=target
                chunk_df['vor_region']=-1
                chunk_df['vor_point']=-1
                # set weight to 1 
                #chunk_df['totalWeight']=1
                chunk_df['totalWeight']=chunk_arr['EventWeight']
                chunk_df['error']=error
                #print(chunk_df)
                data=data.append(chunk_df, ignore_index=True)
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

    # plot the finite line segments
    for simplex in ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex>=0):
            plt.plot(vertices[simplex, 0], vertices[simplex, 1], '-', color='orange', linewidth=2, alpha =0.6)

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
    
    return plt

def eff_error(a, b, err_a, err_b, threshold):
    '''
        effective error assuming 0 correlations 
    '''
    if (b+a) < threshold:
        err = 0
    else:
        err = (a*a*err_a*err_a + b*b*err_b*err_b)/(b+a)
    return err

def get_significance(s, b1, b2, err_b1, err_b2, threshold):
    '''
        asimov significance
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
    
    

