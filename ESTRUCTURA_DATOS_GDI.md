# Documentación y Diccionario de Datos: Pipeline `GDi/` y Miniscope CA1

Este documento detalla la jerarquía completa de archivos, scripts y estructuras de datos utilizando sus **rutas absolutas completas en el sistema** (`/mnt/NAS/...`), cubriendo desde el directorio raíz del almacenamiento experimental hasta los archivos binarios de bajo nivel de análisis de calcio.

Ruta raíz base:
[`/mnt/NAS/Miniscopes/Reg_CA1/GDi`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi)

---

## 1. Nivel Raíz: `/mnt/NAS/Miniscopes/Reg_CA1/GDi`

Contiene la colección de registros con miniscopios inalámbricos/alámbricos (UCLA V4 Miniscope) y cámaras cenitales de comportamiento en la región CA1 del hipocampo, junto con las herramientas centrales de procesamiento cinemático, control de calidad morfométrico y tracking por DeepLabCut (DLC).

### Carpetas de Sujetos / Animales
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R001`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R001): Directorio de cohorte del animal R001.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R002`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R002): Directorio de cohorte del animal R002.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R003`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R003): Directorio de cohorte del animal R003.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26): Cohorte 2026, Ratón 4, Miniscope 4.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26): Cohorte 2026, Ratón 5, Miniscope 5.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26): Cohorte 2026, Ratón 6, Miniscope 6.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R007`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R007): Directorio de cohorte del animal R007.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R008`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R008): Directorio de cohorte del animal R008.

### Scripts de Procesamiento y Control de Calidad en `GDi/`
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales.py):
  Orquestador principal de tracking conductual. Rastrea recursivamente las carpetas de los animales en busca de videos en `My_WebCam/`, verifica si existen inferencias previas de DeepLabCut (`.h5`, `.csv`) y lanza la ejecución sobre los videos pendientes.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales_shuffle2.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales_shuffle2.py):
  Variante de `analizar_animales.py` configurada para ejecutar inferencias con el modelo DeepLabCut entrenado en el shuffle 2.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales_cambiar suffle 1.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/analizar_animales_cambiar%20suffle%201.py):
  Variante de `analizar_animales.py` configurada para alternar la inferencia con el modelo DeepLabCut shuffle 1.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/analisis_geometria.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/analisis_geometria.py):
  Módulo de validación morfométrica que mide la estabilidad anatómica del animal a lo largo de los cuadros del video.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/score_geometria.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/score_geometria.py):
  Cálculo cuantitativo de distancias euclidianas entre pares anatómicos adyacentes (`Nose-Neck`, `Neck-Spine`, `Spine-TailB`, `TailB-TailE`) para computar un *score de distorsión geométrica* que detecta saltos espaciales espurios de DLC.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_v1.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_v1.py):
  Pipeline inicial de filtrado de trayectorias cinemáticas (filtro de mediana e interpolación de coordenadas de baja probabilidad).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_v2_diagnostic.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_v2_diagnostic.py):
  Versión avanzada de filtrado temporal de tracking con generación automática de paneles de diagnóstico y métricas de error.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_mark.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/filter_tracking_mark.py):
  Filtrado y estabilización de las marcas espaciales y objetos colocados dentro de la arena.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/mark_v2_diagnostic.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/mark_v2_diagnostic.py):
  Diagnóstico de estabilidad y deriva espacial de las marcas de los objetos en la arena.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/prueba_filtro_temporal.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/prueba_filtro_temporal.py):
  Script de prueba unitaria y ajuste de parámetros de suavizado temporal sobre series de tiempo de posición $(x, y)$.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/clasificar_eventos_v2.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/clasificar_eventos_v2.py):
  Clasificador de conducta animal basado en cinemática (cálculo de velocidad instantánea, aceleración, orientación de la cabeza y aproximación/interacción con objetos).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/inspeccionar_eventos_v2.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/inspeccionar_eventos_v2.py):
  Herramienta diagnóstica interactiva para inspeccionar visualmente los eventos conductuales detectados.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/visualizar_marcas.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/visualizar_marcas.py):
  Script de renderizado para reproducir el video cenital superponiendo el esqueleto anatómico predicho por DeepLabCut, las marcas filtradas y los polígonos de los objetos.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/Using`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/Using):
  Archivo marcador temporal de control de ejecución (0 bytes).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/__pycache__`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/__pycache__):
  Caché de bytecode compilado de Python (`.pyc`).

---

## 2. Nivel Animal: `/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26`, `R005_M5_26`, `R006_M6_26`

Representan los datos longitudinales completos de cada animal experimental.

### A. Animal R004: [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26)
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_all_days_minian.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_all_days_minian.py): Script de automatización para ejecutar el pipeline de MiniAn secuencialmente en todas las sesiones de R004.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_all_days_minian.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_all_days_minian.log): Archivo de log con la salida y métricas de procesamiento de MiniAn de todos los días.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_minian_batch.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_minian_batch.py): Pipeline de análisis de calcio por lotes con MiniAn (corrección de movimiento, extracción de ROIs y cálculo de fluorescencia).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_cross_registration.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/run_cross_registration.py): Algoritmo de alineación espacial de centroides y registro cruzado de identidad neuronal inter-bloques e inter-días.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/cross_registration_log.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/cross_registration_log.log): Archivo de registro del proceso de cross-registration de R004.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/06_export_mappings_and_C_for_matlab.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/06_export_mappings_and_C_for_matlab.py): Exportador de matrices de fluorescencia $C$ y mapeos celulares a estructuras `.mat` de MATLAB.

### B. Animal R005: [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26)
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian.py): Automatización de MiniAn versión 1 para R005.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian_v2.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian_v2.py): Automatización de MiniAn versión 2 con manejo robusto de excepciones para R005.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_all_days_minian.log): Log consolidado de las corridas de MiniAn en R005.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_minian_batch.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_minian_batch.py): Pipeline por lotes de MiniAn configurado para las grabaciones del miniscopio M5.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_cross_registration.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/run_cross_registration.py): Script de registro celular cruzado entre sesiones de R005.

### C. Animal R006: [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26)
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_all_days_minian.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_all_days_minian.py): Automatizador general de MiniAn para R006.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_all_days_minian.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_all_days_minian.log): Log de ejecución general de MiniAn en R006.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_minian_batch.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_minian_batch.py): Pipeline por lotes de MiniAn para el sensor M6.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/minian_batch.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/minian_batch.log): Log detallado de la ejecución por lotes de MiniAn.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_cross_registration.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_cross_registration.py): Versión 1 del pipeline de registro celular cruzado.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_cross_registration_v2.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/run_cross_registration_v2.py): Versión 2 optimizada del algoritmo de registro celular cruzado.

---

## 3. Nivel Sesión: Rutas de Carpetas Experimentales

### A. Sesiones en `R004_M4_26`
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab C1_2026_06_16`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab%20C1_2026_06_16): Habituación en campo abierto 1.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab C2_2026_06_17`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab%20C2_2026_06_17): Habituación en campo abierto 2.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab L_2026_06_18_stripes_in_11_and_12`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/Hab%20L_2026_06_18_stripes_in_11_and_12): Habituación en pista lineal con rayas visuales en 11 y 12.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19): Test 1 (Same Day, CNO / DREADDs).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T2_SD_VEH__2026_06_21`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T2_SD_VEH__2026_06_21): Test 2 (Same Day, Vehículo control).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T3_XsS_CNO_2026_06_24`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T3_XsS_CNO_2026_06_24): Test 3 (Cross Session, CNO / DREADDs).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T4_XsS_VEH_2026_06_30_repeticion`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T4_XsS_VEH_2026_06_30_repeticion): Test 4 (Cross Session, Vehículo control repetición).

### B. Sesiones en `R005_M5_26`
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab C1_2026_06_16`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab%20C1_2026_06_16): Habituación en campo abierto 1.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab C2_2026_06_17`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab%20C2_2026_06_17): Habituación en campo abierto 2.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab L_2026_06_18`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/Hab%20L_2026_06_18): Habituación en pista lineal.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19): Test 1 (Same Day, Vehículo control).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T2_SD_CNO_2026_06_21`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T2_SD_CNO_2026_06_21): Test 2 (Same Day, CNO / DREADDs).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T3_XsS_VEH_2026_06_24`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T3_XsS_VEH_2026_06_24): Test 3 (Cross Session, Vehículo control).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T4_XsS_CNO_2026_06_24`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T4_XsS_CNO_2026_06_24): Test 4 (Cross Session, CNO / DREADDs).

### C. Sesiones en `R006_M6_26`
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab C1_2026_06_16`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab%20C1_2026_06_16): Habituación en campo abierto 1.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab C2_2026_06_17`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab%20C2_2026_06_17): Habituación en campo abierto 2.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab L_2026_06_18`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/Hab%20L_2026_06_18): Habituación en pista lineal.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T1_SD_VEH_2026_06_20`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T1_SD_VEH_2026_06_20): Test 1 (Same Day, Vehículo control).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T2_SD_CNO_2026_06_22`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T2_SD_CNO_2026_06_22): Test 2 (Same Day, CNO / DREADDs).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T3_XsS_VEH_2026_06_25`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T3_XsS_VEH_2026_06_25): Test 3 (Cross Session, Vehículo control).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T4_XsS_CNO_2026_06_27`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R006_M6_26/T4_XsS_CNO_2026_06_27): Test 4 (Cross Session, CNO / DREADDs).

---

## 4. Archivos a Nivel de Cada Sesión

*(Tomando como ruta representativa: `/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19`)*

- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cents.pkl`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cents.pkl):
  Coordenadas espaciales $(X, Y)$ de los centroides somáticos neuronales detectados en la jornada.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/mappings.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/mappings.csv) y [`mappings.pkl`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/mappings.pkl):
  Matriz de correspondencia y alineación de neuronas a lo largo de los bloques (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/mappings_matlab.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/mappings_matlab.csv):
  Versión con indexación 1-based para su importación directa en MATLAB.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/shiftds.nc`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/shiftds.nc):
  Campo vectorial de desplazamientos espaciales (NetCDF / xarray) calculado para corregir la deriva inter-bloque del miniscopio.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_alignment.png`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_alignment.png):
  Diagnóstico visual que superpone las huellas neuronales alineadas de cada fase con código de colores.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cell_contours.csv.gz`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cell_contours.csv.gz):
  Polígonos vectoriales comprimidos con las coordenadas que delimitan el contorno de cada ROI neuronal.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_contours.pdf`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_contours.pdf):
  Reporte gráfico vectorial en PDF con la distribución espacial de todos los contornos celulares sobre el sensor.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/unit_id_to_C_column.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/unit_id_to_C_column.csv):
  Diccionario que asocia el identificador unificado de la neurona con el índice de columna en la matriz de actividad $C$.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/open_field_alignment.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/open_field_alignment.mat):
  Parámetros de calibración geométrica y límites espaciales de la arena para análisis cinemático en MATLAB.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_minian_batch.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_minian_batch.py):
  Script de ejecución local del pipeline MiniAn para esta sesión.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_cross_registration.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_cross_registration.py):
  Script local de registro cruzado.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/minian_batch.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/minian_batch.log):
  Registro de salida de la corrida de MiniAn de la sesión.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_log.log`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/cross_registration_log.log):
  Registro de salida del cross-registration de la sesión.
- Copias históricas fechadas para reproducibilidad:
  - [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_minian_batch_USED_20260708_212242.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_minian_batch_USED_20260708_212242.py)
  - [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_cross_registration_USED_20260709_023631.py`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/run_cross_registration_USED_20260709_023631.py)

---

## 5. Nivel Bloque: Subsesiones y Datos Raw

*(Tomando como ejemplo el bloque `B_SAMPLE`: `/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE`)*

### Subcarpetas de Fases de la Sesión:
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/A_AV1`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/A_AV1): Ambiente Vacío 1 (registro basal).
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE): Fase de Muestra / Exploración de objetos.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/C_TEST`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/C_TEST): Fase de Prueba / Test de discriminación o memoria.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/D_AV2`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/D_AV2): Ambiente Vacío 2 (registro posterior).

### Contenido Interno de `My_V4_Miniscope/`:
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/0.avi`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/0.avi) (y `1.avi`, `2.avi`, ...): Archivos de video crudo capturados por el sensor CMOS del Miniscope.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/timeStamps.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/timeStamps.csv): Marcas de tiempo de hardware (milisegundos) de cada fotograma de calcio.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/headOrientation.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/headOrientation.csv): Registro inercial IMU tridimensional (roll, pitch, yaw, cuaterniones) de la cabeza del ratón.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/minian`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/minian): Almacenamiento Zarr con arrays intermedios de preprocesamiento, corrección de movimiento y CNMF.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/minian_C_for_matlab.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/minian_C_for_matlab.mat): Matriz de actividad de calcio temporal desconvuelta ($C$) exportada para MATLAB.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/metaData.json`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_V4_Miniscope/metaData.json): Configuración del hardware de miniscopio (potencia de LED, gain, framerate).

### Contenido Interno de `My_WebCam/`:
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0.avi`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0.avi): Video cenital de alta resolución del comportamiento del animal en la arena.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120.h5`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120.h5): Archivo HDF5 de DeepLabCut con las coordenadas $(X, Y)$ y valores de probabilidad de cada parte anatómica y objeto frame a frame.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120.csv): Versión en formato CSV plano de las predicciones de DeepLabCut.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120_meta.pickle`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/My_WebCam/0DLC_HrnetW32_DLC_objectsMay13shuffle2_detector_190_snapshot_120_meta.pickle): Metadatos del modelo DLC, nombres de bodyparts y configuración de inferencia.

### Metadatos a Nivel de Bloque:
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/metaData.json`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/metaData.json): Metadatos generales de la sesión de adquisición de UCLA.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/notes.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/notes.csv): Notas de texto registradas durante el ensayo.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/behaviorTracker`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/behaviorTracker): Parámetros de seguimiento integrados en el hardware.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/experiment`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19/B_SAMPLE/experiment): Parámetros del protocolo de registro de UCLA.

---

## 6. Nivel Especializado: `/mnt/NAS/Miniscopes/Reg_CA1/GDi/<ANIMAL>/<SESION>/CaliAli_results`

La suite **CaliAli** (*Calcium Alignment*) en MATLAB procesa la alineación no rígida entre bloques, la extracción por CNMF y la caracterización morfo-fisiológica de cada neurona individual.

*(Ruta de ejemplo real: `/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results`)*

### 1. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_traces_R005_M5_26_T1_SD_VEH_2026_06_19.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_traces_R005_M5_26_T1_SD_VEH_2026_06_19.mat)
Archivo en formato HDF5 (MATLAB v7.3) que almacena las series de tiempo completas:
- **`traces.C`**: Matriz de trazas de fluorescencia desconvueltas y filtradas (señal de calcio pura sin ruido neuropilar ni fotoblanqueo).
- **`traces.C_raw`**: Matriz de trazas de fluorescencia crudas extraídas directamente de los píxeles de cada ROI.
- **`traces.S`**: Raster de eventos o picos discretos de calcio inferidos (*deconvolved spike events*).
- **`traces.session`**: Vector discreto que identifica a qué sub-bloque (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`) pertenece cada fotograma temporal.
- **`traces.raw_frame`**: Índices originales de los fotogramas en los archivos de video crudo.
- **`traces.stripe_frames`**: Vector de fotogramas marcados por artefactos de ruido de franjas ópticas (*stripes*).
- **`fps`**: Frecuencia de muestreo del registro de fluorescencia en Hertz (típicamente ~30 Hz).
- **`traces_note`**: Documentación técnica del algoritmo sobre las transformaciones aplicadas.

### 2. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.csv`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.csv)
Tabla cuantitativa donde cada fila corresponde a una neurona unificada y alineada:
- **`neuron`**: Identificador entero único de la célula registrada.
- **`centroid_x`, `centroid_y`**: Coordenadas espaciales $(X, Y)$ del centro del soma celular en el plano del sensor.
- **`area_px`**: Superficie del soma de la célula en píxeles.
- **`n_sessions_active`**: Cantidad de sub-bloques (entre 1 y 4) en los que la neurona presentó disparos verificados.
- **Columnas por fase de la tarea** (calculadas para `A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`):
  - `*_events`: Conteo total de eventos de calcio registrados en esa fase.
  - `*_events_min`: Tasa media de disparo (eventos por minuto).
  - `*_snr`: Relación Señal-Ruido (*Signal-to-Noise Ratio*) calculada para esa célula en esa fase.
  - `*_active`: Flag binario (1 o 0) que define si la célula estuvo activa según el umbral mínimo de eventos.

### 3. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.mat)
Estructura en MATLAB que complementa al CSV anterior incorporando:
- **`summary.contours_xy`**: Polígonos vectoriales completos con los vértices perimetrales de cada soma.
- **`summary.contour_level`**: Umbral matemático utilizado para la segmentación del contorno somático.
- **`summary.correlation_image`**: Matriz 2D de la correlación espacial local de intensidades del tejido.
- **`summary.dims_hw`**: Dimensiones del plano focal (alto y ancho del sensor CMOS en píxeles).

### 4. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_options_R005_M5_26_T1_SD_VEH_2026_06_19.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_options_R005_M5_26_T1_SD_VEH_2026_06_19.mat)
Estructura `CaliAli_options` con los hiperparámetros de cómputo:
- **`downsampling`**: Factores de submuestreo espacial y temporal aplicados para optimización computacional.
- **`preprocessing`**: Parámetros de filtrado espacial pasa-banda y sustracción morfológica de fondo.
- **`motion_correction`**: Configuración del algoritmo NoRMCorre para la corrección intra-bloque.
- **`inter_session_alignment`**: Configuración de algoritmos de optimización de traslación y registro afín entre bloques.
- **`cnmf`**: Hiperparámetros de factorización de matrices no negativas (radio neuronal $gSig$, componentes iniciales y umbrales de aceptación).

### 5. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_result_R005_M5_26_T1_SD_VEH_2026_06_19.mat`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_result_R005_M5_26_T1_SD_VEH_2026_06_19.mat)
Archivo binario maestro de gran tamaño (típicamente 200 MB - 850 MB) con el objeto estructurado `neuron`:
- Matrices de huella espacial $A$ (distribución de intensidad en el sensor de cada célula identificada).
- Matrices de fondo espacial $b$ y temporal $f$ (neuropilo circundante difuso).
- Residuos de la deconvolución y métricas de convergencia matemática del modelo CNMF.

### 6. [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19.txt`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19.txt)
Reporte cuantitativo en texto plano con las métricas de registro inter-bloque:
- **`BV score`** (*Background Variation score*): Métrica de estabilidad residual del fondo óptico entre bloques.
- **Tabla comparativa de alineación**: Compara 3 estados:
  1. *Non-aligned*: Sin corrección inter-bloque.
  2. *Translations*: Corrección de desplazamientos rígidos en $X$ e $Y$.
  3. *CaliAli*: Alineamiento elástico no rígido completo.
- **Métricas de salida**:
  - `Correlation scores`: Correlación espacial entre las proyecciones alineadas de cada bloque.
  - `Mean Corr. Score`: Correlación promedio global del campo celular.
  - `Crispness`: Nitidez/contraste espacial de la imagen promedio alineada (a mayor valor, menor desenfoque y mejor superposición).
  *(En sesiones de un solo bloque como `Hab C1` y `Hab C2`, reporta: `"Single session: no inter-session alignment"`).*

### 7. Diagnósticos Visuales (`*.png`) en `CaliAli_results/`
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19_projections.png`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19_projections.png):
  Panel diagnóstico comparativo con las proyecciones compuestas del campo de visión antes y después de la alineación para evaluar rotaciones o derivas.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19_scores.png`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_alignment_R005_M5_26_T1_SD_VEH_2026_06_19_scores.png):
  Curvas de convergencia y evolución de las métricas de correlación espacial y crispness a lo largo de las iteraciones de registro.
- [`/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.png`](file:///mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/T1_SD_VEH_2026_06_19/CaliAli_results/CaliAli_summary_R005_M5_26_T1_SD_VEH_2026_06_19.png):
  Panel principal de control de calidad: integra el mapa de correlación espacial del tejido con los contornos celulares detectados, junto con histogramas de SNR, tasa de eventos y perfiles de actividad a lo largo de las subsesiones.
