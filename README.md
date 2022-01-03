# Bio-Scoring
Desarrollo enfocado en la predicción de compuestos por medio del cáculo de métricas a partir de los datos obtenidos del análisis genómico y de espectrometría. 

## Espectrometría
Empleando como información de entrada el análisis efectuado por [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) por medio su plataforma web, la propuesta enfoca sus esfuerzos en el procesamiento de las predicciones tabuladas en formato html, adjuntas en la carpeta /knownclusterblast/ relativa a cada cepa incluida en el análisis.

  .
  ├── input
  │   ├── antismash
  │   │   ├── CEPA_1
  │   │   │   ├── knownclusterblast 
  │   │   │   │   ├── region1
  │   │   │   │   │   ├── table_1....html
  │   │   │   │   │   ├── table_2....html
  │   │   │   │   │   └── ...
  │   │   │   │   ├── region2
  │   │   │   │   │   ├── table_1....html
  │   │   │   │   │   ├── table_2....html
  │   │   │   │   │   └── ...
  │   │   │   │   └── ...
  │   │   ├── CEPA_2
  │   │   │   ├── knownclusterblast 
  │   │   │   │   ├── region1
  │   │   │   │   │   ├── table_1....html
  │   │   │   │   │   ├── table_2....html
  │   │   │   │   │   └── ...
  │   │   │   │   ├── region2
  │   │   │   │   │   ├── table_1....html
  │   │   │   │   │   ├── table_2....html
  │   │   │   │   │   └── ...
  │   │   │   │   └── ...
  │   │   │   └── ......
  │   │   └── ......
  │   └── ......
  └── ...


## Genómica
Empleando como información de entrada el análisis efectuado por [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) por medio su plataforma web, la propuesta centra sus esfuerzos en la lectura y cálculo de métricas por medio de la comparación con espectrogramas internos. Para evaluar el nivel de correlación entre espectrogramas se empleó la [_Cosine similarity_](https://en.wikipedia.org/wiki/Cosine_similarity). Concretamente, los archivos necesarios poseen como extensión:
  - clustersummary.
  - METABOLOMICS.mgf.

## Flujo de trabajo
Por como esta configurada la propuesta, está efectua la lectura, filtrado y verificación de archivos necesarios dejados en la carpeta /input/ del proyecto, de modo que solo es necesario extraer todos los archivos de [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) en la carpeta /input/ y los de [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) en una subcarpeta /input/antismash/, quedando la siguiente estructura:

  .
  ├── input
  │   ├── antismash
  │   ├── clusterinfo
  │   ├── clusterinfosummary...
  │   │   └── ....clusterinfosummary
  │   ├── gnps_molecular_network_graphml
  │   ├── groupmapping_convertm
  │   ├── qiime2_output
  │   ├── result_specnets_DB
  │   └── METABOLOMICS....mgf
  └── ...

Realizada la lectura de datos, se opera como sigue:

<ol>
  <li>Ordena y agrupa todas las tablas html en un único <i>dataframe</i>, el cual esta constituido por:</li>
  <ul>
    <li>user_BGC</li>
    <li>CEPA</li>
    <li>REGION</li>
    <li>source_bgc_gene</li>
    <li>MIBIG_ID</li>
    <li>identity_percent</li>
    <li>Coverage_percent</li>
    <li>BLAST_Score</li>
    <li>E_value</li>
  </ul>
  
  <li>Se condensa la información genómica a nivel de BGC, promediando las métricas <i>identity_percent</i>, <i>Coverage_percent</i>, <i>BLAST_Score</i> y <i>E_value</i>, añadendo como atributo el número de genes que componen a cada BGC, resultando:
   <ul>
    <li>user_BGC</li>
    <li>CEPA</li>
    <li>REGION</li>
    <li>N_gen</li>
    <li>MIBIG_ID</li>
    <li>identity_percent</li>
    <li>Coverage_percent</li>
    <li>BLAST_Score</li>
    <li>E_value</li>
  </ul>
    
  <li>Se compará el listado de BGC resultantes con los contenidos en <i>matched_mibig_gnps_update.csv</i> (rescatado de <a href="https://github.com/sdrogers/nplinker">NPLinker</a>), filtrando aquellos que no estén contenidos en el archivo interno. Por su parte, los elementos que si están contenidos se les añade los atributos <i>mibig_inchi</i>, <i>mibig_smiles</i> y <i>mibig_name</i> (equivalente al compuesto), resultando en:
  <ul>
    <li>user_BGC</li>
    <li>CEPA</li>
    <li>REGION</li>
    <li>N_gen</li>
    <li>MIBIG_ID</li>
    <li>identity_percent</li>
    <li>Coverage_percent</li>
    <li>BLAST_Score</li>
    <li>E_value</li>
    <li>mibig_inchi</li>
    <li>mibig_smiles</li>
    <li>mibig_name</li>
  </ul>
  </li>
  
  <li>Sujeto a lo anterior, <i>matched_mibig_gnps_update.csv</i> está caracterizado con incluir la relación (link) entre un número acotado de BGC con un número acotado de espectrogramas, indicando el compuesto generado. En consecuencia, aplicado el filtro del punto anterior se establece un grado de similitud entre los BGC del usuario con los de este listado (a través de los atributos <i>identity_percent</i>, <i>Coverage_percent</i>, <i>BLAST_Score</i> y <i>E_value</i>) y en consecuencia con un número acotado de espectrogramas. Bajo esta idea, empleando como métrica <a href="https://en.wikipedia.org/wiki/Cosine_similarity"><i>Cosine similarity</i></a> se establece el nivel de correlación entre los espectrogramas del usuario (provenientes de <a href="https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp">GNPS</a>) con los asociados a los BGC predichos ya filtrados, añadiendo como atributos los asociados al análisis de espectrometría. 
  <ul>
    <li>user_BGC</li>
    <li>CEPA</li>
    <li>REGION</li>
    <li>N_gen</li>
    <li>MIBIG_ID</li>
    <li>identity_percent</li>
    <li>Coverage_percent</li>
    <li>BLAST_Score</li>
    <li>E_value</li>
    <li>mibig_inchi</li>
    <li>mibig_smiles</li>
    <li>mibig_name</li>
    <li>user_clustersumID</li>
    <li>user_SpectraID</li>
    <li>user_Strains</li>
    <li>gnps_key</li>
    <li>gnps_inchi</li>
    <li>gnps_smiles</li>
    <li>gnps_id</li>
    <li>gnps_SpectraID</li>
    <li>simil_cos</li>
    <li>Id</li>
  </ul>
  </li>
  Este <i>dataframe</i> es exportado en formato CSV con el nombre <i>output.csv</i>.
  
    .
    ├── input
    └── output
        └── output.csv
    
  
  <li>Concluido el procesamiento de datos, se procede a establecer el frente de pareto teniendo por defecto la relación <i>identity_percent</i> y <i>simil_cos</i>. El primer término (asociado a los BGC) puede ser modificado al configurar el párametro <i>scoreX</i>, teniendo como opciones:
  <ul>
    <li>Coverage_percent</li>
    <li>BLAST_Score</li>
    <li>E_value</li>
  </ul>
  </li>
  

    
  <li>Concatenado con el punto anterior, se exporta un archivo <i>graphml</i> con el nombre <i>results.graphml</i> </li>
    
    .
    ├── input
    └── output
        ├── output.csv
        └── results.graphml
  
    
    
</ol>

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
    <key attr.name="id" attr.type="string" for="node" id="id"/>
    <key attr.name="product" attr.type="string" for="node" id="product"/>
    <key attr.name="pindex" attr.type="int" for="node" id="pindex"/>
    <key attr.name="label" attr.type="string" for="node" id="label"/>
    <key attr.name="size" attr.type="float" for="node" id="size"/>
    <key attr.name="r" attr.type="int" for="node" id="r"/>
    <key attr.name="g" attr.type="int" for="node" id="g"/>
    <key attr.name="b" attr.type="int" for="node" id="b"/>
    <key attr.name="x" attr.type="float" for="node" id="x"/>
    <key attr.name="y" attr.type="float" for="node" id="y"/>
    <key attr.name="weight" attr.type="double" for="edge" id="weight"/>
    ....
```


## Instalación y puesta en marcha
```python
pip install ...
```

```python
from proposal.graph import CreateGraph 
path  = '..../project/'
# Procesamiento de datos (puntos 1-2-3-4)
graph = CreateGraph(path, MaxCosSimil=True) # default MaxCosSimil=True 
# Exportación de resultados (puntos 5-6)
graph('CEPA', scoreX='identity_percent')    # default scoreX='identity_percent'
```
## Referencias y desarrollos empleados 
La estructura propuesta basa su desarrollo:
- [NPLinker](https://github.com/sdrogers/nplinker) utilizando parte del código y base de datos para el tratamiento de datos de entrada.
- [cdk_pywrapper](https://github.com/sebotic/cdk_pywrapper) utilizando su desarrollo para identificar los _fingerprint_ de los datos rescatados de [NPLinker](https://github.com/sdrogers/nplinker).


