#!/usr/bin/env python
# coding: utf-8

# # Integration postprocessing
# # CODE SNIPPET
# 


# In[1]:


import scanpy as sc
import pandas as pd
from matplotlib import rcParams


# In[2]:


path = ""

# new
adata_normalised_beta = sc.read_h5ad(path + "adata.h5ad")


# In[4]:


# normalize ALL genes (use size factors calculaated on QC genes)
# some of the markers are ambient genes
adata_normalised_beta = adata_normalised_beta.raw.to_adata()
adata_normalised_beta.X /= adata_normalised_beta.obs['size_factors'].values[:,None] # This reshapes the size-factors array
sc.pp.log1p(adata_normalised_beta)


# In[5]:


# gene names and eid gene names
df_genes = pd.read_pickle(path + "gene_names.pkl")


# In[6]:


sc.pl.umap(adata_normalised_beta, color = ['study'], size=10, use_raw=False)


# plotting markers
# In[7]:


markers_selection_disease={
    'healhty': ['XRD2']
 #Also others
}


# In[8]:

def gene_name_to_eid(gene_list):
    print("gene")
    gene_list_eid = [df_genes["GeneEID"][df_genes["GeneName"] == x] for x in gene_list]
    gene_list_eid = [x.tolist()[0] for x in gene_list_eid]
    return gene_list_eid


# In[9]:

new_genes = ['CD81']

markers_selection_new = {
    'genes': new_genes,
}


# In[10]:


for disease in markers_selection_new.keys():
    print(disease)
    genes = markers_selection_new[disease]
    genes_eid = gene_name_to_eid(genes)
    #genes_eid
    print(genes_eid)
    if genes_eid:
        sc.pl.umap(adata_normalised_beta, color = genes_eid, title = genes, size=10)
    else:
        print("No available genes")
        

