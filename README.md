# Bio-Scoring
Desarrollo enfocado en la predicción de compuestos por medio del cáculo de métricas a partir de los datos obtenidos del análisis genómico y de espectrometría. 

## Espectrometría
Empleando como información de entrada el análisis efectuado por [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) por medio su plataforma web, la propuesta enfoca sus esfuerzos en el procesamiento de las predicciones tabuladas en formato html, adjuntas en la carpeta /knownclusterblast/ relativa a cada cepa incluida en el análisis.

## Genómica
Empleando como información de entrada el análisis efectuado por [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) por medio su plataforma web, la propuesta centra sus esfuerzos en la lectura y cálculo de métricas por medio de la comparación con espectrogramas internos de la propuesta. Para evaluar el nivel de correlación entre espectrogramas se empleó la [_Cosine similarity_](https://en.wikipedia.org/wiki/Cosine_similarity). Concretamente, los archivos necesarios poseen como extensión:
  - clustersummary.
  - METABOLOMICS.mgf.


## Flujo de trabajo


## Referencias y desarrollos empleados 
La estructura propuesta basa su desarrollo:
- [NPLinker](https://github.com/sdrogers/nplinker) utilizando parte del código y base de datos para el tratamiento de datos de entrada.
- [cdk_pywrapper](https://github.com/sebotic/cdk_pywrapper) utilizando su desarrollo para identificar los _fingerprint_ de los datos rescatados de [NPLinker](https://github.com/sdrogers/nplinker).


