#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import nbconvert
import pandas as pd
from mprod import table2tensor
from mprod.dimensionality_reduction import TCAM
from sklearn.decomposition import PCA
import seaborn as sn
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import pearsonr
from upsetplot import UpSet
#from skbio.diversity import beta_diversity#Not compatible with current python version 
#from skbio.stats.distance import permanova #Not compatible with current python version 


# In[3]:


def Pseudocount(d):
    p = min(d[d>0]) / 2
    d2 = d + p
    return(d2)

def CLR_norm(x):
    d = np.array(x)
    d = Pseudocount(d)
    d = pd.DataFrame(d.astype(float), columns=list(x.columns))
    step1_1 = d.apply(np.log, axis = 0)
    step1_2 = step1_1.apply(np.average, axis = 1)
    step1_3 = step1_2.apply(np.exp)
    step2 = d.divide(step1_3, axis = 0)
    step3 = step2.apply(np.log, axis = 1)
    return(step3)
def TCA(DF):
    print("Create tensor")
    data_tensor, map1, map3 =  table2tensor(DF, missing_flag=True)
    print("Run TCA")
    tca = TCAM()
    tca_trans = tca.fit_transform(data_tensor)
    print("Generating result tables")
    tca_loadings = tca.mode2_loadings  # Obtain TCAM loadings
    tca_var = tca.explained_variance_ratio_*100 # % explained variation per TCA factor
    tca_df = pd.DataFrame(tca_trans)   # Cast TCA scores to dataframe
    tca_df.rename(index = dict(map(reversed, map1.items())) , inplace = True)    # use the inverse of map1 to denote each row 
                                   # of the TCAM scores with it's subject ID
    return(tca_trans, tca_loadings, tca_var, tca_df)
def Attach_columns(tca_df, Info, Column=["Type"] ):
    Column.append("CS_BABY_BIOME_ID")
    To_add = Info[Column]
    To_add = To_add.drop_duplicates(keep='last')
    tca_df["CS_BABY_BIOME_ID"] = list(tca_df.index)
    tca_df2 =tca_df.merge(To_add, on='CS_BABY_BIOME_ID', how='left')
    return(tca_df2)

def Attach_columns_NEXT(tca_df, Info, Column=["Type"] ):
    Column.append("NEXT_ID")
    To_add = Info[Column]
    To_add = To_add.drop_duplicates(keep='last')
    tca_df["NEXT_ID"] = list(tca_df.index)
    tca_df5 =tca_df.merge(To_add, on='NEXT_ID', how='left')
    return(tca_df5)


def Check_Loads(tca_loadings, Bugs, Top=15, Axis=0, Return=False, Plot=True):
    tca_loadings = pd.DataFrame(tca_loadings, index= Bugs)
    if Plot == False: return(tca_loadings)
    import matplotlib.pyplot as plt
    sn.kdeplot(data=tca_loadings, x=Axis)
    plt.show()
    plt.clf()
    #Biggest contributors to axis
    print(tca_loadings.sort_values(Axis).iloc[0:Top, Axis])
    print("====================================")
    print(tca_loadings.sort_values(Axis, ascending=False).iloc[0:Top,Axis])
    if Return== True:
        return(tca_loadings)

def Count_NA(x):
    x = x[x.notnull()]
    return(x.shape[0])    
def Quantify_Availability(Info, Recoded_time, Keep = [90] ):
    Info["Timepoint2"] =  Recoded_time
    Time_availability = Info.pivot(index='CS_BABY_BIOME_ID', columns='Timepoint_numeric', values='NG_ID')
    #Count total number of non-missing data
    Number_available = Time_availability.apply(Count_NA, axis=0)
    print(Number_available)
    
    #Filter according to keep
    Times_keep = Time_availability[Keep]
    Times_keep = Times_keep.dropna()
    
    Time_availability[ Time_availability.notnull() ] = True
    Time_availability[ Time_availability.isnull() ] = False
    
   
    
    #UpSet Plot
    Time_availability = Time_availability.reset_index()
    Time_availability.set_index( sorted(set(Recoded_time) ) , inplace=True)
    plt = UpSet(Time_availability, subset_size='count',  show_counts=True).plot()
    

    return(Times_keep)

def Cluster(tca_var, tca_df, N_components=2, N_clusters=2 ):
    import matplotlib.pyplot as plt
    print( "Total variability explained: " + str( round(tca_var[0:N_components].sum(), 2) ) + "%" )
    linkage_data = linkage(tca_df.iloc[:,0:N_components] , method='ward', metric='euclidean')
    dendrogram(linkage_data)
    plt.show()
    plt.clf()
    
    hierarchical_cluster = AgglomerativeClustering(n_clusters=N_clusters, affinity='euclidean', linkage='ward')
    labels = hierarchical_cluster.fit_predict(tca_df.iloc[:,0:N_components])
    
    tca_df["Cluster"] = list(labels)
    fig = sn.scatterplot(data=tca_df, x=0, y=1, hue="Cluster")
    fig.set(xlabel= 'Axis 1 ({Perc}%)'.format(Perc=round(tca_var[0],2)), ylabel='Axis 2 ({Perc}%)'.format(Perc=round(tca_var[1],2)) )
    plt.show()
    plt.clf()
    
    return(tca_df)

def Plot_abundances_vs_time(tca_df, Info, DF, Bugs_plot):
    import matplotlib.pyplot as plt
    To_merge = tca_df.reset_index(drop=False)
    To_merge = To_merge[["index", "Cluster"]]
    To_merge = To_merge.rename(columns={"index": "CS_BABY_BIOME_ID"})
    Info_clusters =Info.merge(To_merge, on='CS_BABY_BIOME_ID', how='left')
    Info_clusters = Info_clusters.dropna()
    DF["Cluster"] = list(Info_clusters["Cluster"])
    for Bug in Bugs_plot:
        Bug_abundance = DF[[Bug, "Cluster"]]
        Bug_abundance = Bug_abundance.reset_index(drop=False)
        Bug_short = Bug.split("|s__")[-1]
        fig = sn.boxplot(data=Bug_abundance, y=Bug, x="Timepoint", hue="Cluster", )
        #fig.map(sn.swarmplot, Bug, 'Timepoint')
        #fig.set_ylim( Bug_abundance[Bug].min() , Bug_abundance[Bug].max() )
        fig.set(ylabel=Bug_short + " (Clr)")
        fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()
        plt.clf()

def Do_permanova(Phenos, Pheno_list, tca_df,Perm=2000):
    Results_permanova = {}
    for Pheno in Pheno_list:
        To_add = Phenos[ ["CS_BABY_BIOME_ID", Pheno] ]
        To_add = To_add.dropna()
        To_add = To_add.astype({"CS_BABY_BIOME_ID": str}, errors='raise') 
        To_add = To_add.set_index("CS_BABY_BIOME_ID", drop=True)

        #Remove missing variables
        Not_Missing = []
        for i in list(tca_df.index):
            if i  in list(To_add.index): Not_Missing.append(i)   
        For_distance = tca_df.loc[ Not_Missing ]
        Distance = beta_diversity(metric = 'sqeuclidean', counts =For_distance, ids= For_distance.index)
   
        try:  
            Result = permanova(Distance, To_add, column=Pheno, permutations=Perm)
            P = Result["p-value"]
            Results_permanova[Pheno] = P
        except:
            Results_permanova[Pheno] = np.nan
    return(Results_permanova)


# In[6]:


# Read data, and only keep species 
Info = pd.read_csv("Desktop/CS_Baby_Biome/OLD_ANALYSIS/METADATA_INFANTS_EARLY_CS_BABY_BIOME_09_06_2023_UPDATED_FEEDING.txt", sep="\t")
print (Info.shape)
DF = pd.read_csv("Desktop/CS_Baby_Biome/submission/2024_submission/analysis/species_filtered_early_atleast_2_infants_CS_Baby_Biome_0.05_RA_01_04_2024.txt", sep="\t")
print (DF.shape)
DF['ID'] = DF.index
DF = DF[DF['ID'].isin(Info.NG_ID)]


DF.ID = DF.ID.astype("category")
DF.ID.cat.set_categories(Info.NG_ID, inplace=True)
DF = DF.sort_values(["ID"])

DF2 = DF.drop(["ID"], 1)
DF2 = CLR_norm(DF2) # CLR only on species level 
#DF2 = DF2.drop("UNCLASSIFIED", 1)
print(DF2.shape)


# In[7]:


# format and make into 3D tensor
Info2 = Info[["CS_BABY_BIOME_ID", "Timepoint_numeric"]]
DF3 = pd.concat([Info2.reset_index(drop=True), DF2.reset_index(drop=True)], axis=1)
#First two columns (ID and Timepoint need to become a multi-level index)
Arrays = [ list(DF3.CS_BABY_BIOME_ID), list(DF3.Timepoint_numeric) ]
tuples = list(zip(*Arrays))
index = pd.MultiIndex.from_tuples(tuples, names=["CS_BABY_BIOME_ID", "Timepoint_numeric"])
DF3.index = index
DF3 =DF3.drop(["CS_BABY_BIOME_ID", "Timepoint_numeric"], 1)
print (DF3.shape)


# In[8]:


print(len(DF3.index))


# In[9]:


print (DF3.index)


# In[11]:


tca_trans, tca_loadings, tca_var, tca_df = TCA(DF3)


# In[12]:


tca_df2 = Attach_columns(tca_df, Info, Column=["feeding_mode_pragmatic"] )


# In[13]:


colors = ["olive", "pink", "skyblue"]
fig = sn.scatterplot(data=tca_df2, x=0, y=1, hue="feeding_mode_pragmatic", palette=colors)
fig.set(xlabel= 'Axis 1 ({Perc}%)'.format(Perc=round(tca_var[0],2)), ylabel='Axis 2 ({Perc}%)'.format(Perc=round(tca_var[1],2)) )


# In[14]:


tca_df3 = Attach_columns(tca_df, Info, Column=["Randomization_AB_all"] )


# In[15]:


import seaborn as sn

# create a scatterplot with red and yellow colors
colors = ["gold", "red"]
fig = sn.scatterplot(data=tca_df3, x=0, y=1, hue="Randomization_AB_all", palette=colors)
fig.set(xlabel= 'Axis 1 ({Perc}%)'.format(Perc=round(tca_var[0],2)), ylabel='Axis 2 ({Perc}%)'.format(Perc=round(tca_var[1],2)) )


# In[16]:


tca_df.to_csv("Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TCAM_CS_BABY_BIOME_01_04_2024.txt")


# In[ ]:




