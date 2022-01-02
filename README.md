# Bio-Scoring
Desarrollo enfocado en la predicción de compuestos por medio del cáculo de métricas a partir de los datos obtenidos del análisis genómico y de espectrometría. 

## Espectrometría
Empleando como información de entrada el análisis efectuado por [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) por medio su plataforma web, la propuesta enfoca sus esfuerzos en el procesamiento de las predicciones tabuladas en formato html, adjuntas en la carpeta /knownclusterblast/ relativa a cada cepa incluida en el análisis.

## Genómica
Empleando como información de entrada el análisis efectuado por [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) por medio su plataforma web, la propuesta centra sus esfuerzos en la lectura y cálculo de métricas por medio de la comparación con espectrogramas internos de la propuesta. Para evaluar el nivel de correlación entre espectrogramas se empleó la [_Cosine similarity_](https://en.wikipedia.org/wiki/Cosine_similarity). Concretamente, los archivos necesarios poseen como extensión:
  - clustersummary.
  - METABOLOMICS.mgf.


## Flujo de trabajo
Por como esta configurada la propuesta, está efectua la lectura, filtrado y verificación de archivos necesarios dejados en la carpeta /input/ del proyecto, de modo que solo es necesario extraer todos los archivos de [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) en la carpeta /input/ y los de [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) en una subcarpeta /input/antismash/. Realizada la lectura de datos, la propuesta opera como sigue:

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
    <li></li>
    <li>Id</li>
  </ul>
  </li>
  
  
</ol>



## Referencias y desarrollos empleados 
La estructura propuesta basa su desarrollo:
- [NPLinker](https://github.com/sdrogers/nplinker) utilizando parte del código y base de datos para el tratamiento de datos de entrada.
- [cdk_pywrapper](https://github.com/sebotic/cdk_pywrapper) utilizando su desarrollo para identificar los _fingerprint_ de los datos rescatados de [NPLinker](https://github.com/sdrogers/nplinker).


