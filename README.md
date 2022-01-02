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
  
  <li>Se condensa la información genómica a nivel de BGC, promediando las métricas <i>identity_percent</i>, <i>Coverage_percent</i>, <i>BLAST_Score</i> y <i>E_value</i> y añadendo como atributo el número de genes que componen a cada BGC, resultando:
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
  <li>Milk</li>
</ol>

## Referencias y desarrollos empleados 
La estructura propuesta basa su desarrollo:
- [NPLinker](https://github.com/sdrogers/nplinker) utilizando parte del código y base de datos para el tratamiento de datos de entrada.
- [cdk_pywrapper](https://github.com/sebotic/cdk_pywrapper) utilizando su desarrollo para identificar los _fingerprint_ de los datos rescatados de [NPLinker](https://github.com/sdrogers/nplinker).


