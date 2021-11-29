import pandas as pd
import numpy as np
import os, sys, html5lib
from .parsers.mgf import LoadMGF
from .metabolomics import mols_to_spectra
from .score import fast_cosine, fast_cosine_shift

class processing_metrics:
  def __init__(self, path):
    """Carga de archivos auxiliares y del usuario"""
    #===================================================================================================
    # Carga de archivos auxiliares
    self.internal_map = pd.read_csv(os.getcwd()+'/proposal/data/matched_mibig_gnps_update.csv')
    self.BGCs  = self.internal_map.mibig_id.unique()
    self.GNPSs = self.internal_map['# mgf_spectrum_id'].unique()
    ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([os.getcwd()+'/proposal/data/matched_mibig_gnps_update.mgf'])
    self.spectrum_inter = mols_to_spectra(ms2, metadata)
    print('BGCs en el internal_map:                  ',self.BGCs.shape[0])
    print('Spectrum en el internal_map:              ',self.GNPSs.shape[0])

    #===================================================================================================
    self.clustersum    = None
    self.spectrum_user = None
    self.file_html     = None
    if os.path.isdir(path):
      try:
        # Carga de archivos de espectrometría del usuario
        self.clustersum = [k for k in os.listdir(path) if np.all([os.path.isdir(path+k), 'summary' in k])][0]
        self.clustersum =  pd.read_csv(path+self.clustersum+'/'+os.listdir(path+self.clustersum)[0], sep='\t')
        mgf_file = [k for k in os.listdir(path) if '.mgf' in k][0]
        ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([path+mgf_file])
        self.spectrum_user = mols_to_spectra(ms2, metadata)
      except:
        print("En la ruta {:s} no están contenidos los datos de espectrometria".format(path))
        print("- clustersum")
        print("- mgf_file")
      try:
        # Rutas con datos de BGCs del usuario
        self.file_html = [path+'antismash/'+k+'/knownclusterblast/' for k in os.listdir(path+'antismash')]
      except:
        print("En la ruta {:s} no están contenidos los datos de salida de antiSmash".format(path))
    else:
      print("En la ruta {:s} no existe".format(path))

  def BGC(self, save=True):
    """Generación de datos de BGC a partir de las tablas de salida de antiSmash"""
    _bgc = []
    for path_html in self.file_html:
      FOLDER  = path_html.split('/')[-3]
      folders = [k for k in os.listdir(path_html) if os.path.isdir(path_html+k)]
      files   = [k for k in os.listdir(path_html) if not os.path.isdir(path_html+k)]
      for region in folders:
        for gene in os.listdir(path_html+region):
          source_bgc_gene = gene.split('_mibig')[0]
          txtList=np.array([kx for kx in os.listdir(path_html) if ('.txt' in kx)&(kx.split('_c')[-1].replace('.txt','')==region.split('region')[-1])])
          boolList=[]
          for kk in txtList:
            list_of_lists = []
            with open(path_html+kk) as f:
              list_of_lists.append([source_bgc_gene in line for line in f])
            boolList.append(np.any(list_of_lists))
          user_BGC = txtList[boolList][0].replace('.txt','').replace('.._c','...region00').replace('_c','.region00')
          df = pd.read_html(path_html+region+'/'+gene)[0]
          df = df[df.columns[[0,2,4,5,6,7]]]
          df['source_bgc_gene']=source_bgc_gene
          df['REGION']=region
          df['CEPA']=FOLDER
          df['user_BGC']=user_BGC
          df.rename(columns={'MIBiG Cluster':'MIBIG_ID',
                              '% ID':'identity_percent',
                              '% Coverage':'Coverage_percent',
                              'BLAST Score':'BLAST_Score',
                              'E-value':'E_value'}, inplace=True)
          _bgc.append(df[['user_BGC','CEPA', 'REGION', 'source_bgc_gene','MIBIG_ID','identity_percent','Coverage_percent','BLAST_Score','E_value']])
    df = pd.concat(_bgc).reset_index(drop=True)
    #===================================================================================================
    # Filtrado de BGC predichos por antiSmash con los contenidos en el internal_map
    print('Tamaño original del archivo:              ',df.shape[0])
    print('Total de BGCs en HTML:                    ',df.MIBIG_ID.unique().shape[0])
    userBGC = df[[k in self.BGCs for k in df.MIBIG_ID]]
    print('Total de BGCs contenidos en internal map: ',userBGC.MIBIG_ID.unique().shape[0])
    print('Tamaño del archivo filtrando los BGCs:    ',userBGC.shape[0])
    #===================================================================================================
    # Ponderación de métricas, reducción en la cardinalidad de la información de BGC
    BGCsum = []
    for user_bgc in userBGC.user_BGC.unique():
      for miBIG_bgc in userBGC[userBGC.user_BGC==user_bgc].MIBIG_ID.unique():
        dfx = userBGC[(userBGC.user_BGC==user_bgc)&(userBGC.MIBIG_ID==miBIG_bgc)]
        identity_percent, Coverage_percent, BLAST_Score, E_value = dfx.mean()
        BGCsum.append([user_bgc, np.unique(dfx['CEPA'])[0], np.unique(dfx.REGION)[0], len(dfx), miBIG_bgc, identity_percent, Coverage_percent, BLAST_Score, E_value])
    BGCsum = pd.DataFrame(BGCsum, columns=['user_BGC', 'CEPA', 'REGION', 'N_gen', 'MIBIG_ID', 'identity_percent', 'Coverage_percent', 'BLAST_Score', 'E_value'])
    print('Tamaño del archivo resumen post promedio: ',BGCsum.shape[0])
    return BGCsum

  def OutputFile(self, MaxCosSimil=True, save=True):
    """Concatenación de información y cálculo de métrica de espectrometria"""
    BGCsum = self.BGC(save=False)
    AddProps = []
    for n, bgc in enumerate(BGCsum.MIBIG_ID):
      interBGC = self.internal_map[self.internal_map.mibig_id==bgc][['mibig_name','mibig_inchi','mibig_smiles']].dropna()
      interBGC = interBGC.groupby(['mibig_name','mibig_inchi','mibig_smiles'])
      idx = [list(k) for k in list(interBGC.groups.values())]
      interBGC = interBGC.size().reset_index().rename(columns={0:'count'})
      interBGC['Id']=idx
      for k in interBGC.iterrows():
        # BGC Data
        bgcData  = BGCsum.iloc[n].copy()
        bgcData['mibig_inchi']   = k[1]['mibig_inchi']
        bgcData['mibig_smiles']  = k[1]['mibig_smiles']
        bgcData['mibig_product'] = k[1]['mibig_name']
        bgcData['Id'] = k[1]['Id']
        # Spectrum Data
        interSpectrum=self.internal_map.iloc[k[1]['Id']][['mibig_id','# mgf_spectrum_id','mgf_inchikey','mgf_smiles']]
        CosSimilData = []
        CosSimil=[]
        for n_gnps, GNPSkey in enumerate(interSpectrum['# mgf_spectrum_id']):
          spec_gnps = self.spectrum_inter[np.where(self.GNPSs==GNPSkey)[0][0]]
          for n_user, spec_user in enumerate(self.spectrum_user):
            spec_userData = self.clustersum.iloc[n_user]
            score, matches = fast_cosine(spec_user, spec_gnps, 0.2, 0)
            if MaxCosSimil:
              CosSimil.append(score)
              CosSimilData.append([n_user, spec_user.spectrum_id, spec_userData.UniqueFileSources, interSpectrum.iloc[n_gnps], spec_gnps.id, spec_gnps.spectrum_id, bgcData['Id'][n_gnps]])
            else:
              bgcData2 = bgcData.copy()
              bgcData2['user_clustersumID'] = n_user
              bgcData2['user_SpectraID']    = spec_user.spectrum_id
              bgcData2['user_Strains']      = spec_userData.UniqueFileSources
              bgcData2['gnps_key']          = interSpectrum.iloc[n_gnps]['# mgf_spectrum_id']
              bgcData2['gnps_inchi']        = interSpectrum.iloc[n_gnps]['mgf_inchikey']
              bgcData2['gnps_smiles']       = interSpectrum.iloc[n_gnps]['mgf_smiles']
              bgcData2['gnps_id']           = spec_gnps.id
              bgcData2['gnps_SpectraID']    = spec_gnps.spectrum_id
              bgcData2['simil_cos']         = score
              bgcData2['Id']                = bgcData['Id'][n_gnps]
              AddProps.append(bgcData2)
        if MaxCosSimil:
          argmax = np.argmax(CosSimil)
          bgcData['user_clustersumID'] = CosSimilData[argmax][0]
          bgcData['user_SpectraID']    = CosSimilData[argmax][1]
          bgcData['user_Strains']      = CosSimilData[argmax][2]
          bgcData['gnps_key']          = CosSimilData[argmax][3]['# mgf_spectrum_id']
          bgcData['gnps_inchi']        = CosSimilData[argmax][3]['mgf_inchikey']
          bgcData['gnps_smiles']       = CosSimilData[argmax][3]['mgf_smiles']
          bgcData['gnps_id']           = CosSimilData[argmax][4]
          bgcData['gnps_SpectraID']    = CosSimilData[argmax][5]
          bgcData['simil_cos']         = CosSimil[argmax]
          bgcData['Id']                = CosSimilData[argmax][6]
          AddProps.append(bgcData)
        break
      break
    AddProps = pd.DataFrame(AddProps)
    if save:
      AddProps.to_csv('Output.csv')  
    return AddProps
  
  def CreateGraph(self, MaxCosSimil=True, save=True):
    ## añadir la construcción de grafo acá
    output = self.OutputFile(MaxCosSimil=MaxCosSimil, save=False,)