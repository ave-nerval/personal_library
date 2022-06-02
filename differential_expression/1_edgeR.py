#!/usr/bin/env python
# coding: utf-8

# # DE in beta cells

# In[1]:


# Differential expression (EdgeR)



# In[2]:


import pandas as pd
import scanpy as sc
import numpy as np
import pickle

import matplotlib.pyplot as plt
import seaborn as sb
# import upsetplot as usp
from matplotlib.patches import Patch
# import venn 

from sklearn.preprocessing import minmax_scale,maxabs_scale

# import diffxpy.api as de
# from diffxpy.testing.det import DifferentialExpressionTestWald

from statsmodels.stats.multitest import multipletests

from scipy import sparse

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

from rpy2.robjects.packages import importr
grdevices = importr('grDevices')


# In[3]:


ro.r('library(edgeR)')


# In[4]:


path_data=''
path_save_r=''
path_fig=''


# In[5]:


# Load adata

adata_full=sc.read(path_data+'adata.h5ad')


# In[6]:


# Load gene names
df_genes = pd.read_pickle("gene_names.pkl")

# In[7]:


adata_full.obs.columns


# Use beta cells

# In[8]:

# Ct col
ct_col='cell_type_integrated_subcluster'


# In[9]:


# Subset to endocrine cells
adata = adata_full[adata_full.obs[ct_col]=="beta"]
print(adata.obs[ct_col].unique())


# In[10]:

adata.obs['study_sample_design'] = adata.obs["study_sample"].astype(str) + "_" + adata.obs["design"].astype(str)

# Exclude samples not used in comparisons. Use NOD 8w as healthy adult due to being in same cluster as STZ control.

# In[11]:


adata.obs.design.unique()


# In[12]:


adata.obs.study_sample_design.unique().tolist()



# In[14]:


# Subset data to selected samples
adata = adata[adata.obs['study_sample_design'].isin(groups_map.keys()),:]
print(adata.obs['study_sample_design'].unique())


# In[15]:


adata.shape


# In[16]:


# Add group info (healthy or T2D)
adata.obs['group']=adata.obs.study_sample_design.map(groups_map)


#  #### Creat pseudobulk

# Normalize data for pseudobulk DE

# In[17]:


def get_rawnormalised(adata,sf_col='size_factors',use_log=True,save_nonlog=False, use_raw=True):
    """
    Copy raw data from adata and normalise it.
    :param adata: To normalise. Is copied.
    :param use_raw: If true use adata.raw, else use adata.X
    """
    if use_raw:
        adata_rawnorm = adata.raw.to_adata().copy()
    else:
        adata_rawnorm = adata.copy()
    adata_rawnorm.X /= adata.obs[sf_col].values[:,None] # This reshapes the size-factors array
    if use_log:
        if save_nonlog:
            adata_rawnorm.layers['normalised_counts']=sparse.csr_matrix(np.asarray(
                adata_rawnorm.X.copy()))
        sc.pp.log1p(adata_rawnorm)
    adata_rawnorm.X = sparse.csr_matrix(np.asarray(adata_rawnorm.X))
    return adata_rawnorm


# In[18]:


# Normalise#### Exploring possible covariates
adata_norm=get_rawnormalised(adata.raw.to_adata(),sf_col='size_factors',
                               use_log=False,use_raw=False)


# #### Exploring possible covariates
# 3 diseased and 25 healthy samples

# Age

# In[19]:


adata_norm.obs["age_numeric"] = adata_norm.obs["age"].apply(lambda x: x.split(" ")[0]).astype(int)

df_age = adata_norm.obs[["age_numeric", "group", "file", "study", "sex"]]
df_age = df_age.drop_duplicates(subset = "file")

adata_norm.obs["age_cat"] = pd.cut(adata_norm.obs["age_numeric"], [18, 45, 65], labels = ["young", "old"])

df_age["age_cat"] = pd.cut(df_age["age_numeric"], [18, 45, 65], labels = ["young", "old"])

# age distribution (healthy and diseased)
df_age.pivot(columns="group", values="age_numeric").plot.hist(bins=100)
plt.show()


# In[20]:


get_ipython().run_line_magic('matplotlib', 'inline')
ax = sb.countplot(x = "age_cat", data = df_age, hue = "group")


# In[21]:


adata_norm[adata_norm.obs['group'] == "healthy"].obs.age_numeric.min() # 22
adata_norm[adata_norm.obs['group'] == "healthy"].obs.age_numeric.max() # 65
# healthy data range is between 22 and 65 years


# In[22]:


adata_norm[adata_norm.obs['group'] == "T2D"].obs.age_numeric.min() # 45
adata_norm[adata_norm.obs['group'] == "T2D"].obs.age_numeric.max() # 57
# T2D age only between 45 and 57 years


# Study

# In[23]:


# T2D samples in human_2020 (2) and human_healthy (1)
# healthy samples in 5 studies
ax = sb.countplot(x = "study", data = df_age, hue = "group")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)


# Sex

# In[24]:


ax = sb.countplot(x = "sex", data = df_age, hue = "group")
# in total 17 male and 11 female samples
# healthy: 15 male and 10 female
# T2D: 2 male (human_2020 study) and 1 female (human_healthy study)
# sex is thus confounded within study


# #### Back to pseudobulking

# #### Age as categorical

# In[25]:


# Creat pseudobulk of study-ct-group-sex
xs=[]
obss=[]
for group,data in adata_norm.obs.groupby(['study','group','age_cat', 'study_sample_design']):
    xs.append(np.array(adata_norm[data.index,:].X.sum(axis=0)).ravel())
    # Make obs
    obs={'study':str(group[0]),'group':str(group[1]), 'age': str(group[2])}
    obss.append(obs)
xs=pd.DataFrame(np.array(xs),columns=adata_norm.var_names)
obss=pd.DataFrame(obss)


# In[26]:


obss


# In[27]:


print('xs',xs.shape)
print('obss',obss.shape)


# #### DE

# In[28]:


group='group'


# In[29]:


get_ipython().run_cell_magic('R', '-i xs -i obss -i group', '# Creat object\ny<-DGEList(counts = t(xs),  samples = obss)\nprint(dim(y))')


# In[30]:


get_ipython().run_cell_magic('R', '', '# remove lowly expressed genes\nkeep <- filterByExpr(y, group=y$samples[,group])\ny<-y[keep, , keep.lib.sizes=FALSE]\nprint(dim(y))')


# In[31]:


get_ipython().run_cell_magic('R', '', '# Effective library size\ny <- calcNormFactors(y)')


# In[32]:


condition = "group"

obss[condition].unique()


# In[33]:


# Build design matrix
dmat_loc=pd.DataFrame(index=obss.index)
dmat_loc['Intercept']=1


condition='study'
for val_idx,val in enumerate(
    sorted([cl for cl in obss[condition].unique() if cl!='human_2020'])):
    dmat_loc.loc[obss[condition]==val,condition+'_'+str(val)]=1
    
condition='group'
for val_idx,val in enumerate(
    sorted([cl for cl in obss[condition].unique() if cl!='healthy'])):
    dmat_loc.loc[obss[condition]==val,condition+'_'+str(val)]=1
    
    
condition='age'
for val_idx,val in enumerate(
    sorted([cl for cl in obss[condition].unique() if cl!='young'])):
    dmat_loc.loc[obss[condition]==val,condition+'_'+str(val)]=1    

dmat_loc.fillna(0,inplace=True)
dmat_loc=dmat_loc.astype('float')

print('dmat_loc')
display(dmat_loc)


# In[34]:


# Add design to R and name rows
ro.globalenv['design']=dmat_loc
ro.r('row.names(design)<-row.names(y$samples)')
ro.r('design')


# In[35]:


get_ipython().run_cell_magic('R', '', '# Robust dispersion \ny <- estimateGLMRobustDisp(y, design)')


# In[36]:


get_ipython().run_cell_magic('R', '', '# Fit - GLM not QL as used robust dispersion\nfit <- glmFit(y, design)')


# In[37]:


get_ipython().run_line_magic('R', 'fit$design')


# In[38]:


dmat_loc.columns


# In[39]:


# Get DE tables - compare transitions between clusters
summaries={}
for cl in ['T2D']:
    # Coef2 in python indexing format
    coef2=np.argwhere(dmat_loc.columns=='group_'+str(cl))[0][0]

    print('Cl',cl,'coef2',coef2)
    coef2=coef2+1
    # Coef2 in R indexing format
    print('coef2 R',coef2)
    ro.globalenv['coef2']=coef2
    res=ro.r('glmLRT(fit, coef=coef2)$table')
    # Add padj
    res['padj']=multipletests(res["PValue"].values, alpha=0.05, method='fdr_bh')[1]
    summaries[cl]=res


# #### Save

# In[40]:


file=path_save_r+'beta_study_disease_agecat_'
file


# In[41]:


get_ipython().run_cell_magic('R', '-i file', "# Save fits\nsave(y,fit,file=paste0(file,'fits.RData'))")


# In[42]:


# %%R -i file
# Reload
#load(file=paste0(file,'fits.RData'))


# In[43]:


# Check what to save !!
# Save summary tables
pickle.dump(summaries,open(file+'summaries.pkl','wb'))


# In[44]:


# Reload results
#summaries=pickle.load(open(file+'summaries.pkl','rb'))


# ## Analyse results

# In[45]:


ALFC=1
FDR=0.05


# N DE

# In[46]:


summaries.items()


# In[47]:


# Transitions
de_summary=[]
name = "beta"

n_up=summaries["T2D"].query('padj<@FDR & logFC>@ALFC').shape[0]
n_down=summaries["T2D"].query('padj<@FDR & logFC<-@ALFC').shape[0]
print(name,': fdr',FDR,', abs(lfc) > ',ALFC,':',
  'up',n_up, 'down',n_down)
de_summary.append({'type':name,'direction':'up','n':n_up})
de_summary.append({'type':name,'direction':'down','n':n_down})
de_summary=pd.DataFrame(de_summary)


# In[48]:


de_summary["cell type"] = "beta"


# In[49]:


get_ipython().run_line_magic('matplotlib', 'inline')
sb.catplot(x='cell type',y='n',hue='direction',col='type',data=de_summary,
            kind="bar", height=2.5, aspect=1)


# In[50]:


# set!
'''
genes_shared=set()

genes_res=set(summaries["T2D"].query('padj<0.05 & logFC>2').index)
if len(genes_shared)==0:
    genes_shared=genes_res
else:
    genes_shared=genes_shared&genes_res
'''


# In[51]:


genes_de_up=set(summaries["T2D"].query('padj<0.05 & logFC>@ALFC').index)
genes_de_down = set(summaries["T2D"].query('padj<0.05 & logFC<-@ALFC').index)


# In[52]:


genes_de_up
genes_de_down


# ### Genes expression across cell cls

# In[53]:


adata_rawnorm = adata.raw.to_adata()
adata_rawnorm.X /= adata.obs['size_factors'].values[:,None] # This reshapes the size-factors array

sc.pp.log1p(adata_rawnorm)


# In[54]:


# recalculate clusters on beta cells

sc.pp.neighbors(adata,use_rep='X_integrated')
sc.tl.umap(adata)


# In[55]:


res=0.4
sc.tl.leiden(adata, resolution=res, key_added='leiden_r'+str(res), 
             directed=True, use_weights=True)


# In[56]:


sc.pl.umap(adata,color='leiden_r'+str(res),s=100)


# In[57]:


# Prepare pb for DE genes

# DE genes across all comparisons

cl_col_full = "leiden_r0.4"

# Creat pseudobulk
adata_rawnorm.obs[cl_col_full]=adata.obs[cl_col_full]
xs1=[]
vnames1=[]


# In[58]:


adata_rawnorm.obs[cl_col_full]


# In[59]:


for group,data in adata_rawnorm.obs.groupby(cl_col_full):
    xs1.append(np.array(adata_rawnorm[data.index,:
                                     ].X.mean(axis=0)).ravel())
    # Make obs
    # make sure obss is str not int if clusters
    vnames1.append(group)


# In[60]:


xs1=pd.DataFrame(maxabs_scale(np.array(xs1)),columns=adata_rawnorm.var_names,index=vnames1)
adata_genes=sc.AnnData(xs1.T)
# ALL BETA
# Map of main ct per cl and cl to color based on main ct
ct_cmap = {}
ct_cmap["beta"] = "#b5bbe3"


# In[61]:


diabetes = "T2D"
ct = "beta"

# Heatmap of DE genes expression across cell clusters
for direction,lfc_str in [('up','>'),('down','<-')]:
    name='_'.join([diabetes,ct,direction])
    genes=summaries[diabetes].query('padj<@FDR & logFC'+lfc_str+'@ALFC').index

    # Heatmap
    ad=adata_genes[genes,:]


    x_temp=pd.DataFrame( ad.X.T,
                        index=ad.var_names,columns=ad.obs_names)

    # Cell cl anno
    # Row colors for heatmap, containing cell type info
    row_colors=pd.DataFrame({
                'cell_type':ct_cmap["beta"]
    },index=adata_genes.var_names)

    fg=sb.clustermap(x_temp, 
                  row_colors=row_colors,
                  col_cluster=True,row_cluster=True,
                 xticklabels=False, yticklabels=True)
     # Adds block annotation titles as axis labels
    # legend for cell types
    handles = [Patch(facecolor=c) for c in ct_cmap.values()]
    l2=plt.legend(handles, ct_cmap.keys(),
                  title='Main cell type',
               bbox_to_anchor=(1.25, 0.6), bbox_transform=plt.gcf().transFigure)    
    # legend for endocrine
    handles = [Patch(facecolor=c) for c in ['r','g']]
 
    plt.gca().add_artist(l2)
    fg.fig.suptitle(name)
    del ad
    del x_temp


# In[62]:


# What are those genes??
genes_de_df_up = df_genes[df_genes["GeneEID"].isin(list(genes_de_up))]
genes_de_df_up["direction"] = "up"


# In[63]:


genes_de_df_down = df_genes[df_genes["GeneEID"].isin(list(genes_de_down))]
genes_de_df_down["direction"] = "down"


# In[64]:


genes_de_df = genes_de_df_up.append(genes_de_df_down)


# In[65]:


genes_de_df.shape




