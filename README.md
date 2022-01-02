# Bio-Scoring
Desarrollo enfocado en la predicción de compuestos por medio del cáculo de métricas a partir de los datos obtenidos del análisis genómico y de espectrometría. 

## Espectrometría
Empleando como información de entrada el análisis efectuado por [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) por medio su plataforma web, la propuesta enfoca sus esfuerzos en el procesamiento de las predicciones tabuladas en formato html, adjuntas en la carpeta /knownclusterblast/ relativa a cada cepa incluida en el análisis.

## Genómica
Empleando como información de entrada el análisis efectuado por [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) por medio su plataforma web, la propuesta centra sus esfuerzos en la lectura y cálculo de métricas por medio de la comparación con espectrogramas internos de la propuesta. Para evaluar el nivel de correlación entre espectrogramas se empleó la [_Cosine similarity_](https://en.wikipedia.org/wiki/Cosine_similarity). 



## Flujo de trabajo


