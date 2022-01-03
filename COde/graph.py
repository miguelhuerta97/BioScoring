import scipy.io
import os, sys, html5lib
from .export import processing_metrics
import time as time
import numpy as np
import matplotlib.pyplot as plt

class CreateGraph:
  def __init__(self, ProjectPath, MaxCosSimil=True):
    self.metrics = processing_metrics(ProjectPath=ProjectPath)
    if 'df' in list(vars(self.metrics).keys()):
      print("El archivo 'output.csv' ya esta contenido en la carpeta {:s}".format(self.metrics.OutputPath))
      self.df = self.metrics.df
    else:
      self.df  = self.metrics.OutputFile(MaxCosSimil=MaxCosSimil)
    print('Los strains dispinibles son: ', list(self.df['CEPA'].unique()))
    self.mat = scipy.io.loadmat('./COde/data/matrix_internal_map.mat'.replace('/',os.path.sep))['num'][1:,1:]

  def __call__(self, strain, scoreX='identity_percent'):
    self.BST_Pareto(strain, scoreX=scoreX)
    self.BST_Export()
  

  def BST_Pareto(self, strain, flagProblem=1,  scoreX='identity_percent'):
    if flagProblem==1:
      if scoreX in ['identity_percent', 'Coverage_percent', 'BLAST_Score', 'E_value']:
        t0 = time.time()
        self.dfx = self.df[self.df['CEPA'] == strain]  # Filtrado de datos
        n = self.dfx.shape[0]
        print(' (2) Determinación del frente de Pareto ...')
        print('       nbr. strain data: %d'%(n))
        self.paretoindex = np.ones([n])*np.nan
        if n!=0:
          self.paretovalx = self.dfx[scoreX].values
          self.paretovaly = self.dfx['simil_cos'].values*100
          for i in range(n):
              j = np.setdiff1d(np.arange(0,n), i) 
              valx = self.paretovalx[i]
              valy = self.paretovaly[i]
              self.paretoindex[i] = sum((self.paretovalx[j]>valx)&(self.paretovaly[j]>valy))
              
          print('       nbr. nodes Pareto front: %d'%(sum(self.paretoindex==0)))
          print('       tiempo proc.: %.2fs'%(time.time()-t0))
        else:
          print('ERROR.BST : La cepa (%s) no se encuentra'%(strain))
          print('Probar con: ', list(self.df['CEPA'].unique()))
      else:
        print('ERROR.BST : Flag score X (%s) desconocido !!!'%(scoreX))
    else:
      print('ERROR.BST : Flag análisis (%d) desconocido !!!'%(flagProblem))

  
  def BST_Export(self, paretobias=1, edgebias=0.4, gephiConfig=[100, 50, 1]):
    gephiposscale, gephinodesize, gephiedgescale = gephiConfig
    NAME = self.dfx['user_BGC'].values  
    PRODUCT = self.dfx['mibig_product'].values
    print(' (3) Impresión archivo graphml ...')  
    gmlfn=self.metrics.OutputPath+os.path.sep+'results.graphml'
    print('       escribiendo archivo "%s"\n'%(gmlfn))
    fid = open(gmlfn, "a")
    fid.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    fid.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
    fid.write('    <key attr.name="id" attr.type="string" for="node" id="id"/>\n')
    fid.write('    <key attr.name="product" attr.type="string" for="node" id="product"/>\n')
    fid.write('    <key attr.name="pindex" attr.type="int" for="node" id="pindex"/>\n')
    fid.write('    <key attr.name="label" attr.type="string" for="node" id="label"/>\n')
    fid.write('    <key attr.name="size" attr.type="float" for="node" id="size"/>\n')
    fid.write('    <key attr.name="r" attr.type="int" for="node" id="r"/>\n')
    fid.write('    <key attr.name="g" attr.type="int" for="node" id="g"/>\n')
    fid.write('    <key attr.name="b" attr.type="int" for="node" id="b"/>\n')
    fid.write('    <key attr.name="x" attr.type="float" for="node" id="x"/>\n')
    fid.write('    <key attr.name="y" attr.type="float" for="node" id="y"/>\n')
    fid.write('    <key attr.name="weight" attr.type="double" for="edge" id="weight"/>\n')
    fid.write('    <graph edgedefault="undirected">\n')
    # nodes
    colors = plt.cm.get_cmap("jet", paretobias+1)
    colors = np.array([colors(j) for j in range(paretobias+1)])[:,:-1]*255*2
    idx2plot = np.where(self.paretoindex<=paretobias)[0]+1
    maxpi    = np.max(self.paretoindex)
    for i in idx2plot:
      name=NAME[i-1]
      product=PRODUCT[i-1]
      pidx=self.paretoindex[i-1]
      valx=gephiposscale*self.paretovalx[i-1]
      valy=gephiposscale*self.paretovaly[i-1]
      fid.write('        <node id="%s">\n'%(name))
      fid.write('            <data key="product">%s</data>\n'%(product))
      fid.write('            <data key="pindex">%d</data>\n'%(maxpi-pidx))
      fid.write('            <data key="label">%s</data>\n'%(product))
      fid.write('            <data key="size">%.4f</data>\n'%(gephinodesize))
      fid.write('            <data key="r">%d</data>\n'%(colors[int(pidx),0]))
      fid.write('            <data key="g">%d</data>\n'%(colors[int(pidx),1]))
      fid.write('            <data key="b">%d</data>\n'%(colors[int(pidx),2]))
      fid.write('            <data key="x">%.4f</data>\n'%(valx))
      fid.write('            <data key="y">%.4f</data>\n'%(valy))
      fid.write('        </node>\n')
    #edges
    
    for i in idx2plot:   ###### revisar que coincida con el id real
      for j in idx2plot: ###### revisar que coincida con el id real
        nameI=NAME[i-1]
        nameJ=NAME[j-1]
        we = gephiedgescale*self.mat[self.dfx.iloc[i-1].Id, self.dfx.iloc[j-1].Id]
        if (j!=i)&(we>=edgebias):
          fid.write('        <edge source="%s" target="%s">\n'%(nameI,nameJ))
          fid.write('            <data key="weight">%.4f</data>\n'%(we))
          fid.write('        </edge>\n')  
    fid.write('    </graph>\n')
    fid.write('</graphml>\n')
    fid.close()









