# Formato de los archivos .mat del dataset

## 1. Estructura del dataset

La carpeta principal contiene animales:

- R004
- R005
- R006

Dentro de cada animal hay archivos por fecha, por ejemplo:

- R004/2026_06_16_merged.mat
- R004/2026_06_17_merged.mat
- R004/2026_06_18_merged.mat

Hay dos tipos de archivos:

1. Sesion simple
   - tiene tracking espacial (x, y)
   - tiene actividad neuronal (C, S, t)
   - ejemplo: R004/2026_06_17_merged.mat

2. Archivo mergeado con varias sub-sesiones
   - tambien tiene x, y y actividad
   - pero `pos.x`, `pos.y`, `act.C`, `act.S`, `act.t` aparecen como arrays por sub-sesion
   - ademas tiene `act.mapping`
   - ejemplo: R004/2026_06_18_merged.mat

Tengo que preguntar sobre este file (R004/2026_06_16_merged.mat), no entiendo por que no tiene x e y.

---

## 2. sesion simple

ejemplo:

- R004/2026_06_17_merged.mat

estructura:

- `act`
- `pos`
- `sess`

acceso en python:

```python
import scipy.io
import numpy as np

p = r"R004/2026_06_17_merged.mat"
d = scipy.io.loadmat(p, squeeze_me=True, struct_as_record=False)

x = np.asarray(d["pos"].x)
y = np.asarray(d["pos"].y)
S = np.asarray(d["act"].S)
C = np.asarray(d["act"].C)
t = np.asarray(d["act"].t)
```

### Variables

#### `pos`

```python
x = d["pos"].x      # forma: (n_frames,)
y = d["pos"].y      # forma: (n_frames,)
```


- `vx`, `vy` velocidad
- `v` velocidad total
- `hd` head direction

ejemplo: en R004/2026_06_17_merged.mat,

- `pos.x.shape == (8477,)`
- `pos.y.shape == (8477,)`

la posicion del animal se guarda como un vector por frame.

#### `act`

```python
C = d["act"].C      # forma: (n_frames, n_neuronas)
S = d["act"].S      # forma: (n_frames, n_neuronas)
t = d["act"].t      # forma: (n_frames,)
```

en el ejemplo:

- `act.C.shape == (8477, 510)`
- `act.S.shape == (8477, 510)`
- `act.t.shape == (8477,)`

interpretacion:

- cada fila = un frame
- cada columna = una neurona de esa sesion
- `S[:, j]` = actividad de la neurona j a lo largo del tiempo
- `C[:, j]` = señal continua de fluorescencia
- `t` = tiempo en segundos

la neurona 0 es la primera columna de la matriz:

```python
neurona_0 = S[:, 0]
```

esto es local a esa sesion. no significa que sea la misma neurona que en otra sesion.

#### `sess`

contiene metadata de la sesion:

- `day_name`
- `animal_path`
- `day_path`
- `sess_paths`
- `names`
- `time`

ejemplo:

- `sess.day_name = 'HabC2'`
- `sess.sess_paths = '/.../A_H1'`

esto me dice de que condicion/archivo viene la sesion.

---

## 3. Formato con varias sub-sesiones (merged)

ejemplo:

- R004/2026_06_18_merged.mat

estructura:

- `act`
- `pos`
- `sess`

pero `act` y `pos` no son escalares unicos: son arrays de sub-sesiones.

ejemplo:

```python
x = np.asarray(d["pos"].x)
y = np.asarray(d["pos"].y)
S = np.asarray(d["act"].S)
C = np.asarray(d["act"].C)
t = np.asarray(d["act"].t)
```

estas variables tienen forma `(n_subsesiones,)`, no `(n_frames,)`.

por ejemplo:

- `pos.x.shape == (4,)`
- `pos.y.shape == (4,)`
- `act.S.shape == (4,)`
- `act.C.shape == (4,)`
- `act.t.shape == (4,)`

cada elemento de esas listas es una sesion diferente:

```python
x0 = np.asarray(d["pos"].x[0])
y0 = np.asarray(d["pos"].y[0])
S0 = np.asarray(d["act"].S[0])
```

entonces:

- `x[0]` = posicion de la sub-sesion A
- `x[1]` = posicion de la sub-sesion B
- `x[2]` = posicion de la sub-sesion C
- `x[3]` = posicion de la sub-sesion D

lo mismo para `act.S` y `act.C`:

- `S[0]` es una matriz `(n_frames_A, n_neuronas_A)`
- `S[1]` es una matriz `(n_frames_B, n_neuronas_B)`
- etc.

### `act.mapping`

este campo es el que permite relacionar neuronas entre las sub-sesiones.

ejemplo:

```python
M = np.asarray(d["act"].mapping)
print(M.shape)
print(M[:10])
```

resultado en R004/2026_06_18_merged.mat:

- `act.mapping.shape == (812, 4)`

interpretacion:

- cada fila = una celula global del merged file
- cada columna = una sub-sesion (A, B, C, D)
- el valor en `M[i, j]` es el indice local de esa celula en la sub-sesion j
- `nan` significa que la celula no existe en esa sub-sesion

ejemplo:

```python
[355. 335.  nan 406.]
```

significa:

- la misma celula aparece como 355 en A
- como 335 en B
- no aparece en C
- como 406 en D

esto es la referencia para encontrar neuronas comunes entre sub-sesiones.

#### `sess.names`

en este tipo de archivos, sirve para nombrar las sub-sesiones:

```python
sess.names = ['A', 'B', 'C', 'D']
```

y `sess.sess_paths` guarda las rutas de cada una.

---

## 4. Formato sin posicion. encontre que un file no tiene datos de posicion


- R004/2026_06_16_merged.mat

actividad sin tracking espacial.

```python
import scipy.io
import numpy as np

d = scipy.io.loadmat("R004/2026_06_16_merged.mat", squeeze_me=True, struct_as_record=False)
print(d.keys())

C = np.asarray(d["C"])
S = np.asarray(d["S"])
```

donde:

- `C.shape == (n_frames, n_neuronas)`
- `S.shape == (n_frames, n_neuronas)`
- no hay `pos`
- no hay `sess` ni `act`

este archivo no sirve para mapas espaciales
