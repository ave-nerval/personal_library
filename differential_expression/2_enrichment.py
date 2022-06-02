#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc


import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

import helper as h

import pickle


# In[2]:


ro.r('library("hypeR")')
ro.r("source(paste(Sys.getenv('WSC'),'','helper_hypeR.R',sep=''))")


# In[18]:


path_save_r=''
path_data=''


# In[19]:


adata_full=sc.read_h5ad(path_data+'adata.h5ad')
ct_col='cell_type_integrated_subcluster'
adata = adata_full[adata_full.obs[ct_col]=="beta"]


# In[40]:


genes_df = pickle.load(open(path_save_r + ''+'_de_genes.pkl','rb'))


# In[10]:


genes_df


# In[21]:


# Ref genes
ref=adata.var['GeneID-human_2020'].dropna().tolist()
ro.globalenv['ref']=ref

# Get gene sets
print('MSIGdb version:',ro.r(f'msigdb_version()'))
gene_sets_go=ro.r(f"msigdb_gsets_custom(species='Homo sapiens',category='C5',subcategories=c('GO:BP','GO:CC','GO:MF'),size_range=c(5,500),filter_gene_sets=NULL,background=ref)")
gene_sets_kegg=ro.r(f"msigdb_gsets_custom(species='Homo sapiens',category='C2',subcategories=c('KEGG'),size_range=c(5,500),filter_gene_sets=NULL,background=ref)")
gene_sets_reactome=ro.r(f"msigdb_gsets_custom(species='Homo sapiens',category='C2',subcategories=c('REACTOME'),size_range=c(5,500),filter_gene_sets=NULL,background=ref)")
get_ipython().run_line_magic('R', '-i gene_sets_go -i gene_sets_kegg -i gene_sets_reactome -o gene_sets gene_sets=c(gene_sets_go,gene_sets_kegg,gene_sets_reactome)')
print('N gene sets:',len(gene_sets))
ro.globalenv['gene_sets']=gene_sets


# In[41]:


enrich_fdr=0.25
ro.globalenv['gs_fdr']=enrich_fdr
enrich_datas={}

genes=genes_df['GeneName'].dropna().tolist()


# In[26]:


gene_sets


# In[42]:



# Calculate enrichment
ro.globalenv['genes']=genes
res=ro.r(f'hypeR(signature=as.vector(unlist(genes)),genesets=gene_sets,test = "hypergeometric",background =  as.vector(unlist(ref)),pval = 1,fdr = gs_fdr,plotting = FALSE,quiet = TRUE)')
ro.globalenv['res']=res
data=ro.r(f'res$data')
enrich_datas["human"]=data
print('N enriched gene sets:',data.shape[0])


# In[43]:


if data.shape[0]>0:
    # Plot top enriched gene sets
    print('Top enriched gene sets')
    data['recall']=data['overlap']/data['geneset']
    data['query_size']=len(genes)
    h.plot_enrich(data=data.rename(
        {'label':'name','fdr':'p_value','overlap':'intersection_size'},axis=1),
        n_terms=20, save=False,min_pval=10**-30, max_pval=enrich_fdr,percent_size=True,
           recall_lim=(0,1))
    h.plot_enr_heatmap(data=data,n_gs=None,xticklabels=False,yticklabels=True)
        





