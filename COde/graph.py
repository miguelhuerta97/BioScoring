import scipy.io
import os, sys, html5lib
from .export import processing_metrics
import time as time
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

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

  def __call__(self, strain, scoreX='identity_percent', paretobias=1, edgebias=0.4, gephiConfig=[100, 50, 1],
                    FR=3.5, node_color="jet", edge_color="Reds", figsize=(14,14), show=True):
    self.BST_Pareto(strain, scoreX=scoreX)
    self.BST_Export(paretobias=paretobias, edgebias=edgebias, gephiConfig=gephiConfig)
    if show: self.configGraph(FR=FR, node_color=node_color, edge_color=edge_color, figsize=figsize, gephiConfig=gephiConfig)   

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
    NAME    = self.dfx['user_BGC'].values  
    PRODUCT = self.df['mibig_product'].values # El profe en su doc lo tiene así

    print(' (3) Impresión archivo graphml ...')  
    try: os.remove(self.metrics.OutputPath+os.path.sep+'results.graphml')
    except: pass
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
    print(len(idx2plot))

    for i in idx2plot:
      for j in idx2plot:
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

  def configGraph(self, FR=3.5, node_color="jet", edge_color="Reds", gephiConfig=[100, 50, 5], figsize=(14,14)):
    gephiposscale, gephinodesize, gephiedgescale = gephiConfig
    G=nx.readwrite.graphml.read_graphml(self.metrics.OutputPath+os.path.sep+'results.graphml')

    # Color de etiqueta --> Ranking -> grado (cambiar color a tonos azules)
    mapping = {k[0]: k[1]['label'] for k in G.nodes(data=True)}
    G = nx.relabel_nodes(G, mapping)

    # 1) Aplicar distribución "Fruchterman Reingold" default, escalar y centrar.
    pos = nx.fruchterman_reingold_layout(G,FR)
    for node, (x,y) in pos.items():
      G.nodes[node]['x'] = float(x)
      G.nodes[node]['y'] = float(y)

    # 2) Aplicar apariencia nodos: Color --> Particion -> Product
    product = np.array([k[1]['product'] for k in list(G.nodes(data=True))])
    colors = plt.cm.get_cmap(node_color, len(np.unique(product)))
    colors = np.array([colors(j) for j in range(len(np.unique(product)))])[:,:-1]
    colors = {k:colors[n] for n,k in enumerate(np.unique(product))}
    node_color = np.array([colors[k[1]['product']] for k in list(G.nodes(data=True))])
    
    # Tamaño --> Ranking -> pindex (normalización min-max)
    plist = np.array([k[1]['pindex'] for k in list(G.nodes(data=True))])
    nlist = len(np.unique(plist))
    pmax  = np.nanmax(plist)
    pmin  = np.nanmin(plist)
    sizes = ((plist-pmin)/(pmax-pmin)+1)*gephinodesize

    # 3) Aplicar apariencia aristas: Color --> Ranking -> peso (cambiar color a tonos rojos)
    wlist = np.sort(np.array([k[2]['weight'] for k in list(list(G.edges(data=True)))]))
    cmap = plt.cm.get_cmap(edge_color, len(np.unique(wlist)))#.reversed() 
    colors = np.array([cmap(j) for j in range(len(np.unique(wlist)))])[:,:-1]
    colors = {k:colors[n] for n,k in enumerate(np.unique(wlist))}
    edge_color = np.array([colors[k[2]['weight']] for k in list(list(G.edges(data=True)))])
    edge_width = np.array([k[2]['weight']*gephiedgescale for k in list(list(G.edges(data=True)))])

    # Show graph
    fig=plt.figure(figsize=figsize)
    gs = plt.GridSpec(1, 2, width_ratios=[20, 1])
    ax=plt.subplot( gs[0,0] )
    nx.draw(G,pos=pos,node_size = sizes,node_color=node_color, edge_color=edge_color, with_labels=True, font_color="#1E90FF", width=edge_width, ax=ax)
    # nx.draw(G,pos=pos,node_color=node_color, edge_color=edge_color, with_labels=True, font_color="#1E90FF", width=edge_width, ax=ax)

    ax=plt.subplot( gs[0,1] )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = np.min(wlist), vmax=np.max(wlist)))
    sm._A = []
    plt.colorbar(sm, cax=ax)
    ax.title.set_text('Edge weight')
    plt.show()
