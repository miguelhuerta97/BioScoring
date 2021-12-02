# Bio-Scoring
Desarrollo enfocado a predecir el compuesto empleando métricas desde genómicas y de espectrometría.
Para ello se emplean las predicciones tabuladas (en formato html) de antiSmash en conjunto al archivo de salida de GNPS. Concretamente son necesarios el siguiente listado archivos:
  - Los archivos adjuntos en la carpeta /knownclusterblast/, resultate de antiSmash.
  - El archivo clustersummary resultate de GNPS.
  - El archivo METABOLOMICS.mgf resultate de GNPS.

Cabe señalar que la estructura propuesta basa su desarrollo en [NPLinker](https://github.com/sdrogers/nplinker), utilizando parte del código para el tratamiento de datos de entrada, mientras que la propuesta esta diseñada para efectuar la lectura, filtrado y verificación de archivos necesarios dejándolos alojados en una carpeta única compuesta por la extracción de los archivos .zip resultantes de antiSmash y GNPS.
