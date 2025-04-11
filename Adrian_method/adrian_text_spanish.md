# .UBAfiuba 

FACULTAD DE INGENIERÍA

## Desarrollo de algoritmos para el estudio topológico de flujos caóticos

Adrián Matías Barreal

Director:
Ing. Pablo D. Roca
proca@fi.uba.ar

Codirectora:
Dra. Gisela D. Charó
gcharo@fi.uba.ar

Asesora:
Dra. Denisse Sciamarella
dsciamarella@fi.uba.ar

---

[comment]: # (Page break)
# Agradecimientos 

A la suerte y a mis padres, por darme las oportunidades que otros no han tenido.
A Betiana, por el afecto y las palabras de motivación.
A mis docentes de FIUBA, por darme las herramientas para entender el mundo.
A mis tutores de tesis, por su gran paciencia y dedicación.

Sin ellos este trabajo no sería lo que logró llegar a ser.

Adrián M. Barreal, Diciembre 2023.

---

[comment]: # (Page break)
# Índice general 

Introducción ..... 5
I Estado del arte ..... 8

1. Topología ..... 9
1.1. Introducción al concepto de homología ..... 9
1.2. Topología algebraica y complejos celulares ..... 12
1.3. Grupos de homología ..... 13
2. Análisis topológico de atractores caóticos ..... 16
2.1. El sistema de Lorenz y los atractores caóticos ..... 16
2.2. El sistema de Rössler ..... 19
2.3. Nubes de puntos y variedades enramadas ..... 20
2.4. Template versus simplex ..... 21
2.5. Nubes de puntos de sistemas dinámicos ..... 24
3. Topología computacional ..... 26
3.1. Homología persistente ..... 27
3.2. El método BraMAH ..... 27
II Diseño de la solución ..... 30
4. Estrategias de segmentación ..... 31
4.1. El algoritmo k-means ..... 31
4.2. El espacio de características de los sistemas dinámicos ..... 33
4.3. El grafo de vecindad entre segmentos ..... 37
5. Construcción del complejo celular ..... 40
5.1. Generación del conjunto de datos ..... 40
5.2. Supuestos sobre los datos de entrada ..... 42
5.3. Segmentación ..... 42
5.4. Construcción del 1-complejo ..... 44

---

[comment]: # (Page break)
5.5. Identificación de las 2-fronteras candidatas ..... 46
5.6. Identificación de las 2-celdas ..... 47
5.7. Orientación de las celdas ..... 48
5.8. Optimización aleatoria ..... 49
6. Construcción del simplex ..... 51
6.1. Extensión de un complejo a un simplex ..... 51
6.2. Criterio de orientación de flujo ..... 52
6.3. Implementación ..... 53
III Implementación y evaluación ..... 55
7. Arquitectura de la solución ..... 56
7.1. Introducción ..... 56
7.2. Estructura de la solución ..... 56
7.2.1. Contexto y dependencias ..... 56
7.2.2. Componentes internos ..... 58
7.3. Tipos de datos e interfaz programática ..... 58
7.3.1. La nube de puntos ..... 59
7.3.2. La segmentación ..... 59
7.3.3. El grafo de flujo ..... 60
7.3.4. El complejo celular ..... 61
7.3.5. El templex ..... 61
7.4. Evaluación y resultados ..... 62
7.4.1. Efecto del optimizador ..... 62
7.4.2. Selección del parámetro $k$ ..... 64
7.4.3. Ejemplos de aplicación ..... 66
7.4.4. Dificultades adicionales con el algoritmo ..... 69
7.4.5. Sugerencias para investigaciones posteriores ..... 70
8. Conclusiones generales ..... 72
Apéndices ..... 76
A. Homologías ..... 77
A.1. Definiciones ..... 77
A.2. Cálculo de grupos de homología ..... 77
A.2.1. Grupos de homología del toro ..... 77
A.2.2. Grupos de homología de la esfera ..... 80
A.3. Complejos simpliciales ..... 81
A.3.1. El complejo de Čech ..... 81
A.3.2. El complejo de Vietoris-Rips ..... 81
B. Complejidad algorítmica ..... 83
B.1. Parámetros de entrada ..... 83
B.2. Preprocesamiento de datos ..... 83

---

[comment]: # (Page break)
B.3. Segmentación ..... 84
B.4. Construcción del 1-esqueleto del complejo ..... 85
B.5. Construcción del complejo celular ..... 86
B.6. Optimización ..... 87
B.7. Obtención del simplex ..... 87
B.8. Procedimiento conjunto ..... 88
C. Interfaz programática ..... 90

---

[comment]: # (Page break)
# Introducción 

Los últimos diez años han sido testigos de que la geometría, la topología y los algoritmos forman una potente mezcla de disciplinas con numerosas aplicaciones dentro y fuera del mundo académico [1]. Por otra parte, la teoría de sistemas dinámicos no lineales muestra que dos sistemas cuyas soluciones tienen el mismo comportamiento cualitativo son topológicamente equivalentes. Fue Henri Poincaré [2], un matemático y filósofo francés, quien en 1895 advirtió el modo fundamental en que las soluciones de un sistema dinámico dependen de su topología. La equivalencia topológica se verifica en un espacio que representa los posibles estados del sistema y que se denomina espacio de estados o espacio de fases. Los estados que el sistema visita durante su evolución toman la forma de una nube de puntos temporalmente correlacionados. Un cambio significativo en el comportamiento de las soluciones de un sistema dinámico se traduce inmediatamente en un cambio en la clase de equivalencia topológica. Por ejemplo, si ese sistema representa el clima, entonces ese cambio topológico puede representar una transformación abrupta del clima [3]. De aquí que se vuelva fundamental diseñar estrategias para calcular las propiedades topológicas de un sistema dinámico.

Existen técnicas para calcular propiedades topológicas a partir nubes de puntos en un espacio cualquiera. El objeto matemático que permite extraer las propiedades topológicas de un espacio cualquiera, sin restricciones dimensionales, se conoce en topología algebraica como complejo celular [4]. Sus propiedades topológicas, expresadas por los llamados grupos de homología, describen los agujeros en la nube de puntos, sin considerar -como ocurre para los sistemas dinámicos en el espacio de fases- que cada punto de la nube es un estado del sistema y que una solución numérica de un sistema dinámico puede leerse como una secuencia de puntos ordenados en el tiempo. Interesa pues, describir no solamente la estructura sobre la cual se encuentran los puntos, sino también las formas no equivalentes de circular por dicha estructura de acuerdo al flujo entre estados. El objeto matemático que habilita esta doble descripción ha sido introducido recientemente bajo el nombre de $\operatorname{templex}[5]$. Es un complejo celular con sentido temporal. Su éxito respecto del complejo celular tradicional radica en la posibilidad de diferenciar sistemas dinámicos que son homológicamente equivalentes (la naturaleza de los agujeros en las nubes de puntos respectivas no cambian), y que sin embargo, corresponden a comportamientos dinámicos con características distintas.

Robert Gilmore, uno de los propulsores de la topología como programa elemental en teoría de sistemas dinámicos [6], afirma que el estudio de la equivalencia topológica permite decidir si dos sistemas dinámicos son equivalentes; y en particular, si un modelo desarrollado a partir de datos es una representación exacta de un sistema físico, químico, biológico, climático, etc. Comparar obser-

---

[comment]: # (Page break)
vaciones con modelos, observaciones entre sí (obtenidas por ejemplo mediante distintas técnicas de medición), o modelos de determinados sistemas entre sí, es fundamental en numerosas disciplinas. En el campo de la climatología, por ejemplo, existen más de veinte modelos que se utilizan para hacer simulaciones y proyecciones futuras. Aunque todos tienen una formulación basada en las ecuaciones que rigen los fluidos, difieren en varios aspectos: la modelización de los aerosoles, la dinámica glaciar u oceánica, entre otros procesos. Se suele hablar de una Babel de modelos. ¿Cómo saber qué modelo elegir? ¿Hasta qué punto podemos decir que un modelo es una buena representación de las observaciones? ¿Qué modelo usar para predecir los cambios abruptos que pueden darse en el clima? El templex es la herramienta matemática con potencial para dar respuesta a estas preguntas [7], pero ante todo es preciso poder calcularlo de manera sistemática.

Esta tesis se propone desarrollar y evaluar técnicas algorítmicas capaces de construir un templex, es decir, en un complejo celular adosado a un grafo dirigido, a partir de una nube de puntos asociada a un flujo de un sistema dinámico. Dicho grafo provee las reglas de conexión entre las celdas de mayor dimensión del complejo celular, en función de las secuencias permitidas por el flujo. Tratándose de un concepto matemático nuevo, no existen antecedentes en la literatura de un algoritmo que, dada una nube de puntos, calcule un templex. Las técnicas de cómputo de complejos celulares que fueron diseñadas específicamente para estudiar sistemas no lineales [8, 9], y que han sido aplicadas a problemas de ingeniería tales como el transporte y la mezcla en fluidos incompresibles [10, 11], serán tomadas como punto de partida de este trabajo. La descripción topológica que se alcanza con este enfoque resulta no obstante incompleta, puesto que no se toma en cuenta sentido del flujo entre las celdas del complejo celular.

El algoritmo de construcción del templex será implementado en Python, un lenguaje de programación ampliamente utilizado en la comunidad científica debido a su simplicidad, expresividad y a la amplia disponibilidad de herramientas para el procesamiento y visualización de datos. El software resultante será provisto como un paquete de código abierto, con una interfaz similar a la que provee networkx (un paquete para el trabajo con grafos) en donde las estructuras matemáticas se modelizan como clases. Del desarrollo de estos algoritmos depende la posibilidad de utilizar el templex en un gran número de aplicaciones en ciencia e ingeniería. Si bien esta tesis se limitará al estudio de sistemas dinámicos tridimensionales, se trata del primer paso hacia una generalización algorítmica de un concepto matemático que, por definición, no tiene restricción dimensional. Las técnicas serán puestas a prueba con nubes de puntos que resultan de integrar numéricamente conjuntos de ecuaciones diferenciales no lineales para sistemas dinámicos paradigmáticos con soluciones caóticas.

El desarrollo de esta tesis se enmarca en las líneas de investigación del Centro de Investigaciones del Mar y la Atmósfera (CIMA), instituto compartido entre CONICET y UBA, bajo el eje temático de "Métodos matemáticos para estudios de tiempo y el clima". En ese sentido, el trabajo interdisciplinario con los científicos del instituto resulta fundamental para generar herramientas informáticas de análisis topológico con uso práctico en estas áreas de investigación.

La tesis está organizada de la siguiente forma:

- Capítulo I Estado del arte

---

[comment]: # (Page break)
- Sección 1 - Topología: En esta sección se presentan los fundamentos de la topología matemática y cómo se diferencia de otras ramas como la geometría.
- Sección 2 - Análisis topológico de atractores caóticos: En esta sección se presentan los conceptos fundamentales sobre sistemas dinámicos no lineales, atractores y caos. También se introducen las dificultades inherentes al estudio de aquellos sistemas que presentan comportamiento caótico.
- Sección 3 - Topología computacional: En esta sección se presentan los conceptos básicos de la topología computacional, la rama de la topología enfocada en el cálculo algorítmico en base a conjuntos de datos concretos.


# - Capítulo II Diseño de la solución 

- Sección 4 - Estrategias de segmentación: En esta sección se describen las estrategias de segmentación de conjuntos de datos desarrolladas durante el transcurso de la tesis. Estas estrategias fueron especialmente diseñadas para trabajar sobre conjuntos de datos producidos muestreando sistemas dinámicos no lineales.
- Sección 5 - Construcción del complejo celular: En esta sección se describe paso por paso el algoritmo producido durante la tesis, que permite convertir una nube de puntos en un objeto topológico denominado complejo celular.
- Sección 6 - Construcción del templex: En esta sección se describe una extensión al algoritmo presentado en la sección 5 que permite convertir el complejo celular en un templex.


## - Capítulo III Implementación y evaluación

- Sección 7 - Arquitectura de la solución: En esta sección se detalla la implementación final de las estrategias y algoritmos descriptos en las secciones previas. Se proveen detalles sobre la estructura de la solución y ejemplos de código que muestran cómo hacer uso de la misma. Adicionalmente se realiza un análisis de la implementación, incluyendo aplicación a conjuntos de datos concretos, así como problemáticas y limitaciones identificadas durante el desarrollo y cómo podrían solventarse en trabajos posteriores.
- Sección 8 - Conclusiones generales: Conclusiones generales y cierre.

---

[comment]: # (Page break)
# Capítulo I 

Estado del arte

---

[comment]: # (Page break)
# Sección 1 

## Topología

A diferencia de la geometría, que se enfoca en la medición de cantidades tales como área, volumen o longitud, la topología se enfoca en cambio en la conectividad, los huecos, las torsiones [12] y otras propiedades topológicas que caracterizan la forma de los denominados espacios topológicos.

Un espacio topológico es un conjunto $\mathcal{X}$ dotado con una noción de "cercanía" en base a la cual es posible definir conceptos como los de continuidad y conectividad. La topología como campo de estudio nos ofrece herramientas para caracterizar los espacios topológicos en términos de sus propiedades.

Un espacio topológico $\mathcal{R}$ es topológicamente equivalente a otro espacio topológico $\mathcal{T}$ si existe un homeomorfismo entre ellos. Un homeomorfismo entre $\mathcal{R}$ y $\mathcal{T}$ es una transformación continua con inversa continua. Los homeomorfismos son transformaciones que preservan la estructura del espacio topológico.

La topología es una herramienta clave en el análisis de sistemas dinámicos no lineales, ya que permite analizar los patrones de comportamiento del sistema a un nivel mayor que el de la partícula individual, como se verá en más detalle en la sección 2 .

### 1.1. Introducción al concepto de homología

Un elemento fundamental de la topología moderna es el concepto de homología, que captura la noción de huecos o vacíos en un espacio topológico. La homología es una herramienta matemática que nos permite comparar espacios topológicos e identificar sus características esenciales en forma sistemática.

Sea $\mathcal{C}$ el ciclo de la Figura 1.1a. Podemos observar que $\mathcal{C}$ puede ser reducido mediante un homeomorfismo - una deformación continua que estira, comprime o dobla, pero siempre sin cortar y sin pegar- a un punto $\mathbf{x}$ sin romper el ciclo; es decir, sin que ningún punto de la curva pase por puntos fuera del continuo que en la Figura 1.1a se representa en azul.

---

[comment]: # (Page break)
![chunk1_img-0.jpeg](chunk1_img-0.jpeg)

Figura 1.1: a) Un ciclo $\mathcal{C}$ sobre una superficie en el plano. Los puntos de la superficie se muestran en azul. Aquí $\mathcal{C}$ puede ser deformado en forma continua hasta convertirse en un punto. b) Un ciclo $\mathcal{C}^{\prime}$ sobre una superficie en el plano. Aquí $\mathcal{C}^{\prime}$ no puede ser deformado en forma continua hasta convertirse en un punto sin que algún punto de la curva abandone los confines del espacio topológico.

Por otro lado, sea $\mathcal{C}^{\prime}$ el ciclo de la Figura 1.1b. $\mathcal{C}^{\prime}$ contiene en su interior un vacío que impide que el ciclo pueda ser deformado continuamente hasta convertirse en un punto. La existencia de un ciclo con esta cualidad en el continuo de la Figura 1.1b nos permite afirmar que en dicho espacio topológico existe un hueco. Dos espacios topológicos que poseen huecos similares son topológicamente equivalentes.

Nótese que, en el caso de espacios topológicos de mayor dimensión, pueden existir también huecos de mayor dimensión. No es lo mismo un hueco superficial como aquellos de la Figura 1.1, que un hueco volumétrico como podría ser el que se produce restando un volumen esférico de radio $r_{1}$ a otro volumen esférico concéntrico de radio $r_{2}>r_{1}$. La homología contempla estas diferencias y las utiliza para caracterizar los espacios topológicos y así distinguirlos unos de otros.

De aquí se desprenden definiciones básicas como la de conjunto simplemente conexo, aquel en el que cualquier ciclo es trivial, es decir que puede ser reducido a un punto mediante deformación continua. Una superficie esférica, por ejemplo, es un conjunto simplemente conexo.

Henri Poincaré (1854-1912) observó que para caracterizar apropiadamente una superficie, era importante contar los ciclos independientes no triviales [2]. Para superficies orientables, llamó a este número el número de Betti. Para calcular este número, Poincaré definió una aritmética de ciclos orientados con una operación de suma. Bajo este sistema, dos ciclos $a$ y $b$ son equivalentes $(a \equiv b)$ si uno puede ser deformado continuamente para convertirse en el otro. Luego, los términos $b$ y $-b$ representan el mismo ciclo cerrado, recorrido en sentidos opuestos. Algebraicamente se cumple que $c+(-c)+b \equiv b$; es decir, $-b$ es el inverso aditivo de $b$. La suma de ciclos es conmutativa y asociativa, y la suma de cualquier par de ciclos es también un ciclo. A la suma $a+b$ se la puede interpretar geométricamente como un ciclo que se construye recorriendo primero el ciclo $a$ en sentido positivo, y luego el ciclo $b$ en sentido positivo, o viceversa (Figura 1.2) [12].

El sistema cuenta también con un ciclo nulo 0 tal que $a+(-a) \equiv 0$. En principio, cualquier ciclo trivial es equivalente a 0 , aunque existen también ciclos cero no triviales [12]. En base a esta aritmética es posible también definir, para una superficie dada $\mathcal{S}$, un conjunto de ciclos $\mathcal{B}=$

---

[comment]: # (Page break)
![chunk1_img-1.jpeg](chunk1_img-1.jpeg)

Figura 1.2: Posible interpretación de la suma de dos ciclos $a+b$.
$\left\{c_{1}, c_{2}, \ldots, c_{n}\right\}$ tal que cualquier otro ciclo $x$ en $\mathcal{S}$ sea equivalente a alguna suma de elementos de $\mathcal{B}$ :

$$
x \equiv \sum_{i=1}^{n} r_{i} c_{i}, \quad r_{i} \in \mathbb{Z}
$$

Para superficies orientables, el tamaño de este conjunto mínimo de ciclos no triviales es el número de Betti (unidimensional) de la superficie. Para el caso de las superficies no orientables, como la banda de Möbius (Figura 1.3), es posible observar ocurrencias como la siguiente: existe algún ciclo $a \not \equiv 0$ tal que $a+a \equiv 0$. En otras palabras, existe algún ciclo $a$ no trivial que recorrido $m>1$ veces sucesivas produce un ciclo $m a \equiv 0$. En este caso, el número de Betti unidimensional de la superficie es igual a la cantidad de ciclos independientes que no cumplen esta propiedad. Luego, dado un ciclo $a$ tal que $m a \equiv 0$, con $m$ mínimo y $m>1$, decimos que $m$ es un coeficiente de torsión de la superficie. Para espacios de mayor dimensión, Poincaré definió números de Betti y coeficientes de torsión en base a variedades también de mayor dimensión. Poincaré demostró que los números de Betti y los coeficientes de torsión son invariantes topológicos; es decir, son invariantes ante homeomorfismos. Al igual que la característica de Euler, los números de Betti y los coeficientes de torsión pueden ser utilizados para categorizar espacios topológicos.
![chunk1_img-2.jpeg](chunk1_img-2.jpeg)

Figura 1.3: La banda de Möbius, una superficie no orientable. En una superficie no orientable el vector normal es ambiguo.

Formalmente, los números de Betti pueden obtenerse evaluando el rango de los denominados grupos de homología del espacio topológico. Una definición más precisa de homología y grupos de homología se presentará en la sección 1.3.

---

[comment]: # (Page break)
# 1.2. Topología algebraica y complejos celulares 

Para dotar a la topología de un mayor rigor, Poincaré adoptó una estrategia de análisis basada en los denominados complejos simpliciales. Un complejo simplicial es una estructura matemática que puede ser interpretada como una malla de simplices. Un símplice es la generalización de un triángulo a $n$ dimensiones: un $n$-símplice es la envolvente convexa de $n+1$ vértices. Un 0 -símplice es un punto, un 1-símplice es un segmento de recta, un 2-símplice es un triángulo, un 3-símplice es un tetraedro, y así sucesivamente. En la Figura 1.4 se muestra un modelo de un toro construido conectando vértices mediante aristas, llenando triángulos y dotando adicionalmente a cada uno de estos elementos con información geométrica para construir una estructura en tres dimensiones espaciales.
![chunk1_img-3.jpeg](chunk1_img-3.jpeg)

Figura 1.4: Modelo de un toro construido conectando simplices.

En el campo de la topología, modelizar un espacio topológico en forma de complejo simplicial reduce la tarea de calcular sus características topológicas a un procedimiento algebraico mecánico. Dada una descripción matemática de la estructura de la Figura 1.4, por ejemplo, es posible calcular los números de Betti y los coeficientes de torsión del toro en forma sistemática [12]. Un avance significativo en este proceso de sistematización fue gracias a la matemática Emmy Noether, quien observó que la homología es un caso particular de una estructura algebraica denominada grupo, y que la forma más adecuada de trabajar con homologías sería mediante los denominados grupos de Betti, o grupos de homología [12]. Cada grupo de homología representa una clase de huecos o vacíos en el espacio topológico.

Una construcción más general que el complejo simplicial es el denominado complejo celular. Un complejo celular representa una colección de celdas o espacios $n$-dimensionales conectados o "pegados" por sus fronteras, que son a su vez espacios de dimensión inferior. En la Figura 1.5 se muestra una representación gráfica de un complejo celular construido por el algoritmo que se detallará en la Sección 5.

En un complejo celular general, una $n$-celda es un espacio compacto tal que su interior es homeomórfico a $\mathbb{R}^{n}$ y cuya frontera es homeomórfica a la $(n-1)$-esfera. El complejo de la Figura 1.5 está compuesto de 2-celdas, cada una homeomórfica en su interior a $\mathbb{R}^{2}$ y con frontera homeomórfica a una circunferencia. En [4] se provee el siguiente conjunto de definiciones:
Definición 1 Una n-celda es un conjunto cuyo interior es homeomórfico al disco n-dimensional

---

[comment]: # (Page break)
![chunk1_img-4.jpeg](chunk1_img-4.jpeg)

Figura 1.5: Un complejo celular que modeliza el atractor de Rössler.
$D^{n}=\left\{\mathbf{x} \in \mathbb{R}^{n}:\|\mathbf{x}\|<1\right\}$ con la propiedad adicional de que su frontera debe estar dividida en un conjunto finito de celdas de dimensión inferior, denominadas caras de la n-celda.
Definición 2 Un complejo $\mathcal{K}$ es un conjunto finito de celdas:

$$
\mathcal{K}=\bigcup\{\sigma: \sigma \text { es una celda }\}
$$

tal que, si $\sigma$ es una celda en $\mathcal{K}$, todas las caras de $\sigma$ son celdas de $\mathcal{K}, y$ dadas dos celdas $\sigma$ y $\tau$, los conjuntos de los puntos interiores (no frontera) de $\sigma$ y $\tau$ son disjuntos. La dimensión de $\mathcal{K}$ es la dimensión de la celda de dimensión más alta que contiene.

La definición 2 requiere que dos celdas $\sigma_{1}$ y $\sigma_{2}$ de un complejo celular $\mathcal{K}$ tengan interiores disjuntos. Es decir, dos celdas en un complejo celular pueden estar conectadas a lo sumo por sus fronteras. En la Figura 1.6 se muestran ejemplos de objetos que, si bien están compuestos por símplices, no son complejos celulares.
![chunk1_img-5.jpeg](chunk1_img-5.jpeg)

Figura 1.6: Ejemplos de estructuras que, si bien están compuestas por símplices, no cumplen los requisitos para ser complejos celulares. Esto se debe a que las intersecciones entre los interiores de los símplices son no nulas. Imagen tomada de C. Kinsey, Topology of Surfaces [4, Figura 4.3].

# 1.3. Grupos de homología 

En un complejo celular es fundamental saber cómo están pegadas las celdas entre sí. También debemos tener en cuenta las direcciones de los bordes pegados. En principio, una 1-celda $b$, por

---

[comment]: # (Page break)
definición homeomórfica al intervalo $[0,1]$, tiene una dirección natural dada por $f:[0,1] \rightarrow b$. Desde el punto inicial $f(0)$, la celda se orienta hacia el punto final $f(1)$. Esta parametrización de la celda le otorga una definición natural a la expresión $+b$, a la que podemos entender como la instrucción de desplazarnos sobre la 1-celda desde el punto inicial al punto final en el sentido de su orientación.

Tal como es posible asignarle orientación a una 1-celda, es también posible asignarle orientación a celdas de mayor dimensión. Una 2-celda puede tener una de dos orientaciones: sentido horario o anti-horario.

Definición 3 Un 2-complejo $\mathcal{K}$ es dirigido si a cada 1-celda se le asigna una dirección (desde el punto inicial al punto final), y cada 2-celda se le asigna una dirección (en sentido horario o anti-horario).

Definición 4 Sea $\mathcal{K}$ un complejo dirigido. Una $k$-cadena $C$ en $\mathcal{K}$ es una suma $C=a_{1} \sigma_{1}+a_{2} \sigma_{2}+$ $\cdots+a_{n} \sigma_{n}$ donde $\sigma_{1}, \sigma_{2}, \ldots, \sigma_{n}$ son $k$-celdas en $\mathcal{K}$ y $a_{1}, a_{2}, \ldots, a_{n}$ son enteros. Definimos $0 \sigma=\emptyset$.
Definición 5 Sea $\mathcal{K}$ un complejo dirigido. Llamamos $C_{k}(\mathcal{K})$ al grupo conformado por todas las $k$-cadenas de $\mathcal{K}$, con $k=0,1, \ldots, \operatorname{dim}(\mathcal{K})$.

Definimos el operador lineal borde, también llamado frontera, como $\partial_{k}: C_{k}(\mathcal{K}) \rightarrow C_{k-1}(\mathcal{K})$ como

$$
\partial_{k}(\sigma)=\sum_{i=0}^{k}(-1)^{i} \partial_{i}(\sigma)
$$

Dada una $k$-cadena $\sigma$, el elemento $\partial_{k}(\sigma)$ se llama borde de $\sigma$. Una $k$-cadena con borde nulo no implica que no tiene frontera, sino que sus $k-1$ celdas se anulan. En este caso, decimos que la $k-1$ cadena se llama $(k-1)$-ciclo.

La frontera de una 0 -celda $P$ es el conjunto nulo $\partial P=\emptyset$. La frontera $\partial(b)$ de una 1-celda dirigida $b$ se define como $Q-P$, donde $Q$ y $P$ son el punto final e inicial del recorrido implicado por la orientación de la celda. La frontera de una 2-celda orientada $\sigma$ es la 1-cadena conformada por las 1-celdas en la frontera de $\sigma$. El signo de la 1-celda $b$ en la 1-cadena $\partial(\sigma)$ es positivo si la dirección de $b$ es consistente con la orientación de $\sigma$ (sentido horario o anti-horario), o negativo en caso contrario. Teniendo en cuenta los signos, la 1-cadena $\partial(\sigma)$ describe un recorrido alrededor de $\sigma$ en sentido compatible con la orientación de $\sigma$ (Figura 1.7).

# Definición 6 Grupo de Homología 

Se define el $k$-ésimo grupo de homología como:

$$
H_{k}=\operatorname{Núcleo}\left(\partial_{k}\right) / \operatorname{Imagen}\left(\partial_{k+1}\right)
$$

El Núcleo $\left(\partial_{k}\right)$ es el grupo de todos los $k$-ciclos, mientras que $\operatorname{Imagen}\left(\partial_{k+1}\right)$ es el grupo de los $k$ ciclos que son borde de alguna $k+1$-cadena. Entonces, $\operatorname{Núcleo}\left(\partial_{k}\right) / \operatorname{Imagen}\left(\partial_{k+1}\right)$ sólo conserva los ciclos que no son borde de ninguna $k+1$-cadena, es decir, sólo se conservan los ciclos que representan 'huecos' o 'agujeros'. Los grupos de homología reúnen toda la información esencial que se tiene sobre los complejos: las componentes conexas (agujeros de dimensión 0), los ciclos (agujeros de dimensión 1), las cavidades que encierra (agujeros de dimensión 2 ) y las hipercavidades (agujeros de dimensión $n)$.

---

[comment]: # (Page break)
![chunk1_img-6.jpeg](chunk1_img-6.jpeg)

Figura 1.7: Ejemplo de una 2-celda $\sigma$ orientada en sentido anti-horario. La frontera de $\sigma$ está compuesta por 1-celdas orientadas $a, b, c$ y $d$, cuyas orientaciones se indican en la figura mediante flechas. Las 1-celdas se encadenan para formar la frontera $\partial(\sigma)=a+b-c-d$.

En el Anexo A. 2 se provee la definición formal de grupo de homologías en manera detallada junto con ejemplos concretos de cómo hacer uso de estas definiciones para calcular en forma algebraica los grupos de homología de distintos espacios topológicos.

---

[comment]: # (Page break)
# Sección 2 

## Análisis topológico de atractores caóticos

Fue en la década de los 60 que el matemático y meteorólogo Edward Lorenz (1918-2008) desarrolló un modelo simplificado de tres ecuaciones diferenciales ordinarias para estudiar la convección atmosférica, y escribió un programa de computadora para simular su evolución [13]. El descubrimiento clave de Lorenz respecto a este sistema se dio en 1961 cuando volvió a ejecutar una simulación previa con condiciones iniciales aparentemente iguales, salvo por lo que por entonces habría pasado por un despreciable error de redondeo. Lorenz observó que la simulación produjo resultados radicalmente distintos a los esperados: en el transcurso de dos meses simulados, el error de redondeo había pasado de ser despreciable a ser tan grande como la señal misma. La experiencia de Lorenz [13] dio origen a lo que eventualmente sería conocido como la teoría del caos, el estudio de los sistemas dinámicos altamente sensibles a las condiciones iniciales. Lo notable de los sistemas caóticos es que el determinismo no implica predictibilidad. Puesto en las palabras de Lorenz:
"Caos: cuando el presente determina el futuro, pero el presente aproximado no determina aproximadamente el futuro."

### 2.1. El sistema de Lorenz y los atractores caóticos

El sistema de ecuaciones de Lorenz parametrizado por $\sigma, \rho, \beta \in \mathbb{R}$ se escribe de la siguiente manera:

$$
\left\{\begin{array}{l}
\frac{\mathrm{d} x}{\mathrm{~d} t}=\sigma(y-x) \\
\frac{\mathrm{d} y}{\mathrm{~d} t}=x(\rho-z)-y \\
\frac{\mathrm{d} z}{\mathrm{~d} t}=x y-\beta z
\end{array}\right.
$$

Dada una condición inicial $\left(x\left(t_{0}\right), y\left(t_{0}\right), z\left(t_{0}\right)\right)$ en un instante de tiempo inicial $t_{0}$ y dados los parámetros $\sigma=10, \rho=28$, y $\beta=8 / 3$, la resolución del sistema de Lorenz produce una sucesión de estados

---

[comment]: # (Page break)
![chunk1_img-7.jpeg](chunk1_img-7.jpeg)

Figura 2.1: Trayectoria en el espacio de fases obtenida integrando numéricamente el sistema de Lorenz con parámetros $\sigma=10, \rho=28, y \beta=8 / 3$ y condiciones iniciales $\mathbf{p}_{0}=(0,1,1.05)$. El punto inicial de la trayectoria se indicó con un marcador.
que describe la evolución del sistema en el tiempo, dando lugar a lo que se conoce con el nombre de trayectoria en el espacio de fases (Figura 2.1). El espacio de fases de un sistema es un espacio matemático en el que se representan todos los estados posibles del mismo. En la Figura 2.2 se muestran dos trayectorias cuya diferencia en las condiciones iniciales es en términos absolutos pequeña, pero se puede observar que esta diferencia se amplifica significativamente con el paso del tiempo.

El experimento de Lorenz nos permite observar una característica desafortunada de los sistemas dinámicos no lineales: en la práctica, el comportamiento a largo plazo de una partícula en régimen caótico es en general impredecible en forma precisa. Esta impredictibilidad es esencial y no salvable: incluso conociendo con exactitud las reglas que rigen el comportamiento del sistema, los pequeños errores propios del proceso de medición -o la resolución finita en una integración numérica- serán eventualmente amplificados y la predicción podrá diferir significativamente de lo observado. Esta propiedad inherente a la sensibilidad a las condiciones iniciales no impide estudiar y entender el comportamiento global de un sistema operando en régimen caótico. Para este análisis se adoptan técnicas basadas en el estudio de los denominados atractores del sistema dinámico.

Coloquialmente, un atractor es un conjunto de estados al que el sistema tiende, y en el que permanece una vez alcanzado. En general, toda trayectoria suficientemente próxima a un atractor permanecerá próxima incluso ante perturbaciones leves. En la Figura 2.3 se muestra una trayectoria que forma el denominado atractor de Lorenz, el conjunto de estados a los que el sistema de Lorenz tiende dados los parámetros y condiciones iniciales adecuadas. La existencia de los atractores y su manifestación en los sistemas físicos es importante para la teoría del caos porque implica la existencia de patrones subyacentes al comportamiento del sistema. Son estos patrones los que nos permiten adquirir un entendimiento del fenómeno físico, y son el foco de estudio de esta teoría.

---

[comment]: # (Page break)
![chunk1_img-8.jpeg](chunk1_img-8.jpeg)

Figura 2.2: Dos trayectorias en el espacio de fases obtenidas simulando la evolución del sistema de Lorenz desde $t=0$ hasta $t=2.5$ con condiciones iniciales $\mathbf{p}_{0 a}=(1,1,20)$ para la curva azul y $\mathbf{p}_{0 b}=(1.05,1.05,20)$ para la curva roja. La diferencia inicial en las condiciones iniciales es en términos absolutos pequeña, pero se amplifica significativamente con el paso del tiempo. Los estados iniciales de cada trayectoria se indican con marcadores circulares (si bien parece un único marcador), y los estados finales con marcadores en forma de cruz.
![chunk1_img-9.jpeg](chunk1_img-9.jpeg)

Figura 2.3: El atractor de Lorenz. La curva es una trayectoria obtenida integrando el sistema de Lorenz en forma numérica. Los parámetros utilizados fueron $\sigma=10, \rho=28, y \beta=8 / 3$, y las condiciones iniciales $\mathbf{p}_{0}=(1,1,1)$.

---

[comment]: # (Page break)
![chunk1_img-10.jpeg](chunk1_img-10.jpeg)

Figura 2.4: Proceso de formación de un atractor caótico de tipo Rössler espiral. En forma cíclica el espacio se estira, se comprime, se dobla sobre sí mismo y se reinserta en el proceso nuevamente. El resultado es un espacio compuesto por infinitas capas apiladas una sobre otra. Las trayectorias tienen estructura fractal.

Existen diversos tipos de atractores; para un sistema dinámico dado pueden existir conjuntos de estados que tienden eventualmente hacia un punto fijo; pueden existir también estados inestables que tienden a la divergencia; alternativamente, un estado inicial puede evolucionar hasta converger eventualmente a una órbita periódica. El atractor de Lorenz no responde a ninguno de estos patrones, sino que pertenece a una familia de atractores denominados atractores extraños. Ninguna trayectoria sobre un atractor extraño es periódica. Una vez alcanzado el atractor, el sistema circulará indefinidamente en forma no periódica, sin diverger ni converger a ningún punto, dentro de la región finita del espacio de fases en la que atractor extraño yace. Este comportamiento se debe a la naturaleza fractal de los atractores extraños. Un fractal es una entidad geométrica creada en base a patrones similares repetidos a distintas escalas. La geometría fractal de los atractores extraños surge debido a la naturaleza de los procesos que los generan. La analogía utilizada típicamente es la de una masa a la que repetidamente se la estira y se la dobla sobre sí misma como si fuese una frazada (Figuras 2.4 y 2.5). El proceso se repite infinitas veces en forma iterativa dentro de los confines del atractor. Es por ello que la Figura 2.3 es en realidad engañosa: en las regiones en las que las trayectorias parecen describir una superficie plana, hay en realidad infinitas capas superficiales apiladas una encima de la otra. Esto es algo que Lorenz ya había observado y explicado en su paper de 1963, a pesar de que por entonces la idea no tuvo mucho impacto.

# 2.2. El sistema de Rössler 

El sistema de Rössler [14] es un sistema no lineal de ecuaciones diferenciales propuesto por Otto Rössler en 1976 como un modelo simplificado del sistema de Lorenz, con objeto de hacer más accesible el entendimiento cualitativo de los sistemas caóticos. Las ecuaciones que describen la

---

[comment]: # (Page break)
![chunk1_img-11.jpeg](chunk1_img-11.jpeg)

Figura 2.5: Trayectoria en el espacio de fases de un sistema de Rössler producida integrando las ecuaciones del sistema desde $t=0$ hasta $t=600$ con parámetros $(a, b, c)=(0.2,0.2,5.7) y$ condiciones iniciales $\mathbf{p}_{0}=(1,1,1)$.
evolución de este sistema son las siguientes:

$$
\left\{\begin{array}{l}
\frac{d x}{d t}=-y-z \\
\frac{d y}{d t}=x+a y \\
\frac{d z}{d t}=b+z(x-c)
\end{array}\right.
$$

donde $a, b$ y $c$ son parámetros reales. La simplicidad relativa de este sistema está en que dos de sus tres ecuaciones son lineales, a diferencia de lo que ocurre con el sistema de Lorenz en el que sólo una de tres ecuaciones es lineal. En la Figura 2.5 se muestra una trayectoria producida integrando el sistema de Rössler. Al atractor de la Figura 2.5 se lo llama atractor de Rössler. El atractor de Rössler es conveniente por ser relativamente simple pero al mismo tiempo presentar varias de las características que hacen a su análisis un problema de interés para esta tesis.

# 2.3. Nubes de puntos y variedades enramadas 

Desde el punto de vista de la teoría del caos, la topología resulta clave para el estudio de los atractores caóticos porque nos brinda información sobre la estructura global del atractor, imprescindible cuando la trayectoria es impredecible en términos prácticos más que a corto plazo. Resultará deseable entonces contar con procedimientos que nos permitan extraer la información topológica subyacente a un atractor dado.

---

[comment]: # (Page break)
Llamamos nube de puntos a un conjunto $\mathcal{X} \subset \mathbb{R}^{n}$ de vectores posición en algún espacio $n$ dimensional. Una nube de puntos, en el contexto de los sistemas dinámicos, puede ser producida mediante un proceso de muestreo de un sistema físico o generada mediante simulación. La nube de puntos de la Figura 2.6 fue producida tomando muestras de una trayectoria en el atractor de Lorenz, generada mediante la resolución numérica del sistema de Lorenz.
![chunk1_img-12.jpeg](chunk1_img-12.jpeg)

Figura 2.6: Una nube de puntos generada muestreando una instancia del sistema de Lorenz.
Dada una representación gráfica de una nube de puntos, como la de la Figura 2.6, el cerebro humano es capaz de inferir intuitivamente la topología del atractor subyacente a la nube, especialmente en las regiones virtualmente planas (en este caso las "alas" del atractor). Existen sin embargo regiones de más difícil interpretación, tanto intuitiva como algorítmica: se trata de aquellos subespacios donde el flujo converge y aparenta intersecarse a sí mismo. Este tipo de espacio topológico es lo que se denomina una variedad enramada, y sus ramas marcan los sitios donde las trayectorias se separan y parten hacia regiones diferenciadas del espacio de fases. El estiramiento y compresión del flujo es lo que permite la separación de las trayectorias. La cantidad y ubicación relativa de las ramas de una variedad enramada dan lugar a distintos tipos de comportamiento caótico, permitiendo distinguir entre distintos tipos de atractores. En las variedades enramadas, identificar correctamente la topología en aquellos sitios donde confluyen las ramas será la principal complicación al momento de construir el complejo en forma algorítmica en secciones posteriores.

# 2.4. Template versus templex 

El templex [15] es una estructura matemática recientemente desarrollada en el contexto del análisis topológico de sistemas dinámicos no lineales en el espacio de fases. Se trata de un objeto que consiste en un complejo celular dotado de información adicional, un grafo de flujo entre celdas, que permite distinguir entre atractores caóticos que no pueden diferenciarse en términos de los grupos de homología.

La palabra templex es una contracción de las palabras template y complex, aunque también podría asociarse al término temporal complex, debido a que se dota al complejo de una información

---

[comment]: # (Page break)
de índole temporal, vinculada al flujo subyacente. Complex hace referencia al complejo celular tal como se introdujo en la sección 1.2. Un template o knot-holder, por otro lado, es un objeto introducido por los matemáticos Joan S. Birman y Robert F. Williams en 1983 [16]. Fue definido con el objeto de describir el modo en que las trayectorias tridimensionales se anudan en un atractor caótico, haciendo uso de la teoría de nudos [17].

En la figura 2.7 se representa gráficamente un template que describe un atractor que encierra un foco. En este atractor el flujo se abre en dos ramas (o strips) 1 y 2 . La rama 2 presenta una torsión, similar a una banda de Möbius. La rama 1 continúa sin torsión. Las ramas reintegran su flujo en una denominada línea de juntura. Podemos visualizar un pequeño entorno alrededor de la línea de juntura como un 2-complejo celular con tres 2-celdas unidas por una misma 1-celda. Los templates pueden ser representados matemáticamente en términos de matrices de enlace y de juntura, que describen cómo las distintas ramas se tuercen, se entrecruzan y se unen.
![chunk1_img-13.jpeg](chunk1_img-13.jpeg)

Figura 2.7: Un template que describe un atractor que ramifica en dos ramas 1 y 2. La rama 2 presenta una torsión. Las ramas confluyen en una línea de juntura.

La mayor limitación del template es que sólo puede ser construido para espacios de fases cuya dimensión no supere 3. Esta es una limitación inherente a la teoría de nudos, ya que los nudos se deshacen en más de tres dimensiones. Los complejos celulares, por otro lado, tienen la ventaja de permitir representar espacios más generales. Usar complejos y grupos de homología para extraer las propiedades topológicas de nubes de puntos en el espacio de fases asociadas a sistemas dinámicos, en más de tres dimensiones, ha sido un procedimiento exitoso en muchas aplicaciones concretas $[8,9,10,11]$.

El complex demuestra superioridad ante el template no solamente por la posibilidad de extender el análisis topológico a espacios de fases de dimensiones más altas, sino también porque su estudio no depende de la reconstrucción de trayectorias a partir de nubes de puntos para formar los nudos. Esto último puede resultar problemático en el contexto de un experimento real debido al ruido presente en la señal, la cual suele adicionalmente ser de duración relativamente corta. A diferencia del template, las celdas del complejo se construyen sobre conjuntos de puntos pertenecientes a distintas trayectorias, sin necesidad de aislar las trayectorias para observar el modo en que se

---

[comment]: # (Page break)
![chunk1_img-14.jpeg](chunk1_img-14.jpeg)

Figura 2.8: Complejo celular orientado que describe la estructura topológica del atractor de Rössler. El complejo está compuesto por 2-celdas $\gamma_{1}, \gamma_{2}, \ldots, \gamma_{6}$, 0-celdas $0,1,2, \ldots, 9$, y 1-celdas generadas conectando 0celdas unidas en el gráfico por una arista. La 1-celda $(0,1)$ es una denominada línea de juntura, caracterizada por ser parte de la frontera de más de dos 2-celdas ( $\gamma_{2}, \gamma_{5}$ y $\gamma_{6}$ en este caso). Las direcciones de las flechas indican cómo deben unirse o "pegarse" las correspondientes instancias de la 1-celda $(0,1)$ para construir el objeto topológico; las flechas deben pegarse tal que la dirección de las flechas sea consistente. Imagen tomada de Charó, Letellier y Sciamarella, 2022. Templex: a bridge between homologies and templates for chaotic attractors [15, Figura 12].
anudan. Más aún, como afirma la topóloga L. Christine Kinsey, hay mucha más información en un complejo celular que la que expresan sus grupos de homología. Las líneas de juntura de un template, por ejemplo, pueden ser obtenidas en base al complejo observando los bordes comunes entre celdas.

La dificultad que surge, no obstante, al querer reemplazar un template por un complex para aproximar una variedad enramada reside en que la relación entre los conceptos de "rama" en un template y de "agujero" en un complex no es directa. Consideremos el ejemplo de la Figura 2.8 en la que se representa un atractor de Rössler mediante un complejo celular. La 1-celda $\langle 0,1\rangle$ es una línea de juntura en la que convergen los flujos de dos ramas definidas por las secuencias de 2-celdas $\gamma_{1}, \gamma_{2}, \gamma_{3}, \gamma_{5}, \gamma_{1}$ por un lado, y $\gamma_{1}, \gamma_{2}, \gamma_{4}, \gamma_{6}, \gamma_{1}$ por otro. Este complejo tiene un único 0 -generador, $\mathcal{H}_{0}=[[<0>]]$, es decir, tiene una única componente conexa. El único 1-generador es

$$
\mathcal{H}_{1}=[[\langle 0,2\rangle-\langle 0,7\rangle+\langle 2,4\rangle+\langle 4,7\rangle]]
$$

que identifica el agujero central del atractor. No hay traza alguna, en las homologías, de la existencia de dos ramas. Esto es un problema porque en términos de sus grupos de homología el atractor de Rössler resulta indistinguible de un flujo cilíndrico. Es decir, el análisis de los grupos de homología del complejo no pone en evidencia la cualidad característica que distingue al sistema de Rössler de sistemas más simples.

Es allí donde surge el complex, una extensión del concepto de complejo celular que incorpora información sobre las junturas y sobre el flujo que conecta las celdas del complejo formando ciclos (circuitos cerrados o semi-cerrados de celdas) desde y hasta las líneas de juntura, si es que existen.

---

[comment]: # (Page break)
![chunk1_img-15.jpeg](chunk1_img-15.jpeg)

Figura 2.9: El digrafo que describe el flujo entre celdas sobre el complejo celular de la Figura 2.8. Imagen tomada de Charó, Letellier y Sciamarella, 2022. Templex: a bridge between homologies and templates for chaotic attractors [15, Figura 12].

Los ciclos no redundantes desde las junturas de un templex serán equivalentes a las ramas (strips) de un template y recibirán el nombre de stripexes. Esta información se incorpora mediante un grafo dirigido, o digrafo, que describe las transiciones válidas (el flujo) entre celdas. El digrafo de la Figura 2.9 describe el flujo en el espacio de fases entre las celdas del complejo de la Figura 2.8. Este digrafo indica la existencia de dos stripexes dados por las siguientes secuencias de 2-celdas del complejo o nodos del digrafo:

$$
\begin{aligned}
& \mathcal{S}_{1}: \gamma_{1} \rightarrow \gamma_{2} \rightarrow \gamma_{4} \rightarrow \gamma_{6} \rightarrow \gamma_{1} \\
& \mathcal{S}_{2}: \gamma_{1} \rightarrow \gamma_{2} \rightarrow \gamma_{3} \rightarrow \gamma_{5} \rightarrow \gamma_{1}
\end{aligned}
$$

una de las cuales $\left(\mathcal{S}_{2}\right)$ sufre medio giro local, llamado local twist.
El digrafo con el complejo celular conforman el antedicho templex y en conjunto contienen la información necesaria para caracterizar con mejorada precisión al sistema dinámico. Un flujo cilíndrico, por ejemplo, puede ser representado mediante un templex sin stripexes, en contraposición a los dos que caracterizan al atractor de Rössler. Toda esta información puede extraerse mediante operaciones algebraicas a partir del complejo y el digrafo asociado a las celdas del complejo. El templex permite así reunir la información espacio-temporal topológica del sistema en espacios de fases de dimensión arbitraria y expresarla en un conjunto de propiedades que no dependerán del número de celdas del templex construido.

# 2.5. Nubes de puntos de sistemas dinámicos 

Los algoritmos producidos para esta tesis trabajan sobre un tipo particular de nube de puntos, más adecuado al contexto de los sistemas dinámicos determinísticos. Estos algoritmos operan sobre nubes de puntos $\mathcal{C}=\left\{\mathbf{x}_{t} \in \mathbb{R}^{n}: 1 \leq t \leq m\right\}$, con $t, m \in \mathbb{N}$ y $\mathbf{x}_{t}$ una sucesión de puntos ordenada temporalmente. Es decir, no estamos lidiando meramente con una nube obtenida muestreando un espacio topológico en forma arbitraria, sino con una nube obtenida muestreando la traza de la secuencia de estados del sistema que recorre (con el correr el tiempo) una solución de un sistema dinámico. La nube de puntos de la Figura 2.6 fue generada muestreando la trayectoria de un sistema de Lorenz: el conjunto de datos cuenta entonces con un ordenamiento temporal implícito en el orden de las muestras en el archivo de entrada. Haciendo uso del ordenamiento es posible inferir propiedades del flujo en el atractor, lo que resulta útil al momento de construir un templex

---

[comment]: # (Page break)
en base a una nube. Esto no es algo común y no se tienen referencias de algoritmo alguno que haga uso de tal información.

---

[comment]: # (Page break)
# Sección 3 

## Topología computacional

Los métodos que estudian las estructuras de grandes conjuntos de datos enfocándose en la conectividad de dicho conjunto, se engloban en el área denominada Análisis Topológico de Datos $[18,19,20]$. El método más popular se denomina homologías persistentes [21, 22] y tiene el objetivo de encontrar estructuras y patrones en los datos. Este método se utiliza cada vez más en diversos ámbitos de investigación, como por ejemplo en el relacionado con el procesamiento de imágenes [23].

Otro método que permite realizar análisis topológico de datos para variedades enramadas es el denominado Branched Manifold Analysis Through Homologies (BraMAH, Sección 3.2) [8, 9], que permite aproximar una variedad enramada $m$-dimensional inmersa en un espacio de dimensión $n, n \geq m$ con un $m$-complejo celular tal que el número de $m$-celdas es mucho menor que el número de puntos en la nube bajo análisis. Las propiedades topológicas de este tipo de complejo han sido empleadas en ciencias de la ingeniería para identificar con éxito las regiones de no mezclado en un fluido a partir de series temporales de partículas aisladas [24, 11].

Una representación discreta y vectorial como lo es la nube de puntos es una entrada adecuada para algoritmos del campo de la topología computacional. La topología computacional es un campo de estudio que combina conceptos de topología, geometría y ciencias de la computación para desarrollar algoritmos y herramientas para el análisis topológico de datos. La calidad de los resultados producidos por estos algoritmos está fuertemente condicionada por la calidad de los datos de entrada: densidad espacial de muestras no uniforme, un muestreo a frecuencia insuficiente, etc.

El foco de esta tesis está puesto en el desarrollo de algoritmos para la construcción de un $\operatorname{templex}$ a partir de nubes de puntos en el espacio de fases. Más precisamente, el objetivo será construir un complejo celular y un grafo dirigido asociado, en base a una nube de puntos que es una representación discreta de un espacio topológico, de forma tal que el flujo subyacente pueda ser analizado en forma sistemática, junto con la topología de la estructura espacial. En esta sección, resumiremos las características de los métodos mencionados anteriormente: homologías persistentes y BraMAH.

---

[comment]: # (Page break)
# 3.1. Homología persistente 

Dada una nube de puntos, el método de homologías persistentes construye una representación denominada código de barras, característica de la topología subyacente a la nube. Este método se basa en la construcción de una sucesión de complejos simpliciales anidados dependientes de un parámetro de escala $\epsilon$. Dos complejos simpliciales típicos son el complejo de Čech y el complejo de Vietoris-Rips [25]. Ambos complejos dependen de un denominado parámetro de escala $\epsilon \in \mathbb{R}$ que permite definir una noción de cercanía entre puntos. En el Anexo A se definen ambos complejos en términos de sus respectivos procedimientos de construcción.

La estrategia de homologías persistentes no es aproximar variedades enramadas subyacentes, sino construir familias de complejos y analizar rasgos persistentes a través de las escalas. El supuesto fundamental es que, dada una familia paramétrica de complejos $\mathcal{K}_{\epsilon}$ anidados, las características topológicas persistentes para rangos relativamente grandes del parámetro $\epsilon$ son más significativas que aquellas con vida relativamente corta, y es más probable que representen verdaderas características topológicas de $\mathcal{T}$. En la Figura 3.1 se observa la familia de complejos anidados $\mathcal{K}_{\epsilon_{1}} \subset \mathcal{K}_{\epsilon_{2}} \cdots \subset \mathcal{K}_{\epsilon_{7}}$ con $\epsilon_{1} \leq \epsilon_{2} \leq \cdots \leq \epsilon_{7}$ y el código de barras. Cada barra horizontal representa alguna característica topológica en alguno de los grupos de homología $H_{0}, H_{1}, H_{2}, \ldots$, y la longitud de la barra representa la persistencia de la característica en cuestión: una barra de mayor longitud representa un parámetro que existe en un rango más grande de $\epsilon$. El gráfico de la Figura 3.1 [25, Figura 4] muestra códigos de barras para $H_{0}, H_{1}$ y $H_{2}$, generados evaluando las características topológicas de complejos generados para distintos valores de $\epsilon$. En general no hay garantía de que exista algún $\epsilon$ tal que $\mathcal{K}_{\epsilon}$ sea topológicamente equivalente al espacio $\mathcal{T}$ subyacente a la nube, es decir, que los rasgos persistentes no permiten seleccionar, entre los complejos calculados, un complejo celular topológicamente fiel a la estructura de la variedad enramada.

El código de barras en homologías persistentes es una herramienta útil para caracterizar la topología $\mathcal{T}$ subyacente a la nube. La limitación más evidente, no obstante, es que no se trata de una descripción precisa de $\mathcal{T}$, sino de una vía indirecta para caracterizar sus propiedades. Si pudiésemos conocer $\mathcal{T}$ con más precisión podríamos comprender mejor el conjunto de datos y, para el caso puntual de los sistemas dinámicos, el comportamiento del sistema que lo generó.

### 3.2. El método BraMAH

Branched Manifold Analysis through Homologies (BraMAH) [8, 9] es un caso concreto de un método que, dada una nube de puntos como la de la figura 2.6, busca construir un complejo representativo de la variedad enramada en la que yacen los puntos, para determinar -a partir de sus propiedades- los invariantes topológicos del espacio subyacente. Una diferencia importante entre BraMAH y el estudio de las homologías persistentes es que el primero supone que la nube de puntos yace en una variedad enramada y se vale de esta hipótesis para generar un único complejo celular del que se extrae la información topológica mediante procedimientos algebraicos.

Se entiende aquí por variedad enramada $m$-dimensional un espacio topológico tal que cada punto tiene una vecindad topológicamente equivalente a un disco abierto $m$-dimensional. Se dice que tal variedad es Hausdorff si y sólo si dos puntos distintos cualesquiera tienen vecindades disjuntas.

---

[comment]: # (Page break)
![chunk1_img-16.jpeg](chunk1_img-16.jpeg)

Figura 3.1: Robert Ghrist 2008, Barcodes: The persistent topology of data [25, Figura 4]. Códigos de barras para los grupos de homología $H_{0}, H_{1}$ y $H_{2}$ (Sección 1.3) que describen la topología de la nube que se muestra en la parte superior de la imagen. Las barras describen cómo distintas características topológicas aparecen $y$ mueren en función de $\epsilon$.

La segunda condición no se cumple precisamente en la unión entre ramas, es decir, en los lugares que describen el estiramiento y el aplastamiento de un flujo en el espacio de fases. Una variedad enramada es, por tanto, una variedad a la que no se exige que cumpla la propiedad de Hausdorff. Si bien no es habitual en topología algebraica considerar espacios de este tipo ([4], página 65), nada impide construir un complejo celular que lo describa.

Sea $\mathcal{C} \in \mathbb{R}^{n}$ una nube de puntos representativa de un atractor que yace en una variedad enramada de dimensión local $m$, los puntos pueden agruparse en los denominados "parches" que van a ir cubriendo progresivamente la variedad. Los parches están asociados cada uno a una $m$-celda que debe aproximar un $m$-disco (un 3 -disco es una bola rellena). Un parche $\mathcal{P}$ se construye de la siguiente manera: sea $\mathcal{P}\left(x_{0}\right)$ un subconjunto de $n_{c}$ puntos de la nube en torno a un punto $x_{0}$. Para verificar si $\mathcal{P}\left(x_{0}\right)$ puede ser aproximado por una $m$-celda, se analiza si los puntos de $\mathcal{P}$ que pertenecen a una bola con centro $\mathbf{x}_{0}$ y radio $r$ aproximan localmente a un disco $m$-dimensional. La distribución de puntos en $\mathcal{P}\left(x_{0}\right)$ puede analizarse construyendo la matriz $X \in \mathbb{R}^{n_{c} \times n}$ :

$$
X_{i, j}=\frac{1}{n_{c}^{1 / 2}}\left(x_{i, j}-x_{0, j}\right)
$$

donde $n_{c}$ es el número de puntos que se pretende agrupar en una sola $m$-celda.

---

[comment]: # (Page break)
Los vectores singulares de $X$ proveen un sistema de coordenadas local centrado en $\mathbf{x}_{0}$, y sus valores singulares describen la distribución de los puntos en torno a $\mathbf{x}_{0}$. Si $\mathcal{P}\left(x_{0}\right)$ aproxima un disco de dimensión $m$ en $\mathbb{R}^{n}$, deben existir $m$ valores singulares de $X$ que son función lineal del número de puntos en $\mathcal{P}\left(x_{0}\right)$ ordenados de manera creciente según la distancia a $x_{0}$ [26]. Los valores singulares restantes $(n-m)$ medirán la desviación del espacio tangente al disco y escalarán como $r^{\ell}$ con $\ell \geq 2$.

Dado un parche $\mathcal{P}_{k}$ construimos una $m$-celda de la siguiente manera: para cada parche $\mathcal{P} \neq \mathcal{P}_{k}$ tal que $\mathcal{P}_{k} \cap \mathcal{P} \neq \emptyset$, se seleccionan $m$ puntos de $\mathcal{P}$. Estos $m$ puntos representan la intersección de $\mathcal{P}_{k}$ con $\mathcal{P}$. A todos $\operatorname{los} \mathbf{x} \in \mathcal{P}_{k}$ representantes de intersecciones con otros parches se los proyecta sobre el plano de ajuste definido por los puntos de $\mathcal{P}_{k}$, y se los ordena mediante un barrido antihorario respecto de su centro de masa. Esta lista de vértices ordenados define la $m$-celda correspondiente a $\mathcal{P}_{k}$. El conjunto de todas las celdas define el complejo celular. Finalmente, contando con el complejo celular, se extraen las características topológicas del complejo usando técnicas algebraicas.

En resumen, para una nube de puntos dada, BraMAH construye un único complejo sobre el cual basa su análisis algebraico. En principio, la técnica permitirá obtener una descripción topológicamente fiel a la variedad enramada subyacente. Sin embargo, y a diferencia de lo que ocurre con las homologías persistentes donde se construye una infinidad de complejos con características topológicas variables, la calidad del complejo construido con BraMAH resultará crítica para garantizar resultados correctos.

La eventual formación de falsos agujeros (huecos artificiales producidos por las limitaciones del algoritmo) constituye una dificultad que acecha en cualquier método de construcción de complejos celulares a partir de una nube de puntos. Las homologías persistentes contrarrestan esta dificultad con la noción de persistencia. BraMAH, en cambio, no tiene una estrategia algorítmica compensatoria. Veremos no obstante, que esta dificultad se verá mitigada al construir un templex, donde la existencia de falsos agujeros o el relleno inadecuado de un agujero con una celda, no alterará la detección de caminos no equivalentes sobre la estructura.

La motivación principal detrás de esta tesis fue desarrollar un algoritmo para generar un templex a partir de una nube de puntos, asumiendo que se busca aproximar no solamente una variedad enramada, como en BraMAH, sino también un flujo sobre la variedad en cuestión. Como se verá más adelante, ambas hipótesis serán fundamentales en la concepción del conjunto de algoritmos desarrollados en esta tesis.

---

[comment]: # (Page break)
# Capítulo II 

## Diseño de la solución

---

[comment]: # (Page break)
# Sección 4 

## Estrategias de segmentación

Para construir un complejo celular que aproxime una variedad enramada localmente $m$-dimensional, la implementación original de BraMAH aplica una técnica de segmentación, utilizando parches o segmentos que se asemejan a un $m$-disco. Durante esta tesis se evaluaron estrategias de segmentación alternativas basadas principalmente en el algoritmo k-means, en busca de un método robusto que permita construir un complejo que pueda ser luego extendido a un simplex.

### 4.1. El algoritmo k-means

El método k-means es un algoritmo de segmentación ampliamente utilizado en ciencia de datos que, dados un parámetro entero $k>0$ y un conjunto de $n \geq k$ muestras $\mathcal{C} \subset \mathbb{R}^{m}$, produce un conjunto de $1 \leq q \leq k$ segmentos o clusters no vacíos $\left\{\mathcal{S}_{i}\right\}$ tal que

$$
\mathcal{C}=\bigcup_{i=1}^{q} \mathcal{S}_{i} \quad \text { y } \quad i \neq j \Rightarrow \mathcal{S}_{i} \cap \mathcal{S}_{j}=\emptyset
$$

A diferencia de la estrategia originalmente adoptada por BRAMAH que busca la superposición de parches para contar con un registro del modo en que las celdas están unidas, k-means asigna cada punto de la nube a un único cluster. El algoritmo k-means funciona de la siguiente forma:

1. Inicialización: Dada la cantidad de clusters $k$, se inicializa aleatoriamente un conjunto de $k$ centroides, uno para cada cluster. Los centroides iniciales son puntos seleccionados aleatoriamente de entre el conjunto de datos. Llamamos $\mathbf{c}_{i} \in \mathbb{R}^{m}$ al centroide correspondiente al cluster $\mathcal{S}_{i}$.
2. Asignación: Se asigna cada punto de $\mathcal{C}$ al cluster correspondiente a su centroide más cercano acorde a alguna métrica de distancia en $\mathbb{R}^{m}$.
3. Actualización: Para cada cluster $\mathcal{S}_{i}$ se recomputa el centroide $c_{i}$ teniendo en cuenta la asignación de puntos generada en el paso anterior.

---

[comment]: # (Page break)
4. Repetición: Se repite el procedimiento desde el paso 2 hasta la convergencia. La convergencia puede determinarse cuando ya no hay cambios en los clusters de un paso al siguiente, o cuando la distancia de $c_{i}$ en el paso $j$ hacia $c_{i}$ en el paso $j+1$ es menor a un cierto $\epsilon$. También es usual limitar la cantidad de iteraciones que el algoritmo ejecuta.

La decisión de adoptar k-means se debió a las siguientes cualidades del algoritmo:

- Se trata de un algoritmo simple, eficiente y popular, para el que existen implementaciones robustas en Python, el lenguaje adoptado para la implementación de los algoritmos producidos para esta tesis.
- Es flexible. La métrica de distancia a utilizar es arbitraria, así como el espacio en el que está embebido la nube.
- k-means guarda una una relación interesante con los diagramas de Voronoi. En el plano, un diagrama de Voronoi describe un complejo celular donde cada celda se asemeja a un 2-disco (Figura 4.1). Esto es esencialmente lo que la implementación original de BraMAH busca. Cuando la nube de puntos puede ser embebida en el plano, la segmentación producida por k-means se asemeja a un diagrama de Voronoi (Figura 4.2).
![chunk1_img-17.jpeg](chunk1_img-17.jpeg)

Figura 4.1: Diagrama de Voronoi en una región del plano.
Formalmente, un diagrama de Voronoi se define de la siguiente forma. Sea $\mathcal{X}$ un espacio métrico con una función de distancia $d$. Sea $\left\{\mathbf{c}_{i}\right\} \subset \mathcal{X}$ un conjunto finito y discreto de elementos de $\mathcal{X}$. La celda de Voronoi $\mathcal{V}_{i}$ asociada al punto $\mathbf{c}_{i} \in \mathcal{V}_{i}$ es el conjunto de puntos en $\mathcal{X}$ tal que, para todo $\mathbf{p} \in \mathcal{V}_{i}$, la distancia $d\left(\mathbf{p}, \mathbf{c}_{i}\right)$ es menor o igual a cualquier otra distancia $d\left(\mathbf{p}, \mathbf{c}_{j}\right), i \neq j$. No es difícil ver porqué los resultados producidos por k-means se asemejan a un diagrama de Voronoi: el resultado producido por k-means asigna a cada punto de la nube $\mathcal{C}$ a su centroide más cercano, produciendo una partición análoga a la del diagrama de Voronoi.

---

[comment]: # (Page break)
![chunk1_img-18.jpeg](chunk1_img-18.jpeg)

Figura 4.2: Una lámina tridimensional segmentada usando k-means. El resultado se asemeja a un diagrama de Voronoi.

# 4.2. El espacio de características de los sistemas dinámicos 

Sea $\mathcal{C} \subset \mathbb{R}^{3}$ una nube de puntos en $\mathbb{R}^{3}$ producida muestreando el estado de un sistema dinámico confinado dentro del atractor de Rössler (Figura 4.3). Cada uno de los puntos de la nube corresponde a un vector posición del sistema en el espacio de estados. Supongamos entonces que usamos k-means para producir una segmentación del atractor, tal como se hizo con la lámina de la Figura 4.2. Un resultado posible es el que se muestra en la Figura 4.3. Se observa inmediatamente que esta segmentación podría no ser buena desde el punto de vista de la topología: el hueco central del atractor, una característica topológica esencial, no se refleja claramente en la segmentación. Consideremos a modo de contraste la segmentación que se muestra en la Figura 4.4; observaremos allí que el hueco central del atractor no queda contenido en el interior de ningún segmento sino que está rodeado de segmentos, lo cual puede eventualmente facilitar su detección. La segmentación de la Figura 4.4 se obtuvo aplicando k-means sobre una representación de la nube en un espacio de características alternativo a aquel de los vectores posición usado para producir la figura 4.3. La idea de espacio de características se describe a continuación.

El concepto de espacio de características es relativamente popular en el mundo de la ciencia de datos y el aprendizaje automático. La idea fundamental es que virtualmente cualquier entidad descriptible en términos de $m$ cantidades cuantificables o categóricas puede ser representado como un vector en algún espacio $m$ dimensional. Los algoritmos como k-means producen distintos resultados dependiendo de las características seleccionadas para representar las entidades con las que se está trabajando. El siguiente fragmento de código muestra cómo se puede aplicar k-means a una nube de puntos usando el lenguaje de programación Python y el paquete Scikit Learn. El fragmento asume que contamos con una matriz features de $n \times m$ elementos donde la $i$-ésima fila es un vector que representa la $i$-ésima muestra de una nube de puntos en un espacio de características

---

[comment]: # (Page break)
![chunk1_img-19.jpeg](chunk1_img-19.jpeg)

Figura 4.3: Nube de Rössler segmentada en el espacio de vectores posición usando k-means. Esta segmentación pierde el hueco central del atractor.
$m$-dimensional:

```
def clusterize(features):
    # Se estandarizan los datos restando la media y dividiendo
    # por el desvío estandar en cada componente.
    scaler = preprocessing. StandardScaler().fit(features)
    scaled = scaler.transform(features)
    # Se aplica k-means con inicialización aleatoria, k clusters
    # iniciales y un máximo de max_iter iteraciones.
    kmeans = clustering.KMeans(init="random", n_clusters=k, max_iter=max_iter)
    kmeans.fit(scaled)
    # Se retorna un vector donde el i-ésimo elemento es un
    # entero que identifica unívocamente al segmento que contiene
    # a la i-ésima muestra.
    return kmeans.predict(scaled)
```

Se puede ver que el resultado de la segmentación depende del espacio de características en el que se decide representar la nube. La segmentación de la Figura 4.4 fue producida representando cada punto de la nube no en términos de su posición espacial, sino en términos de la dirección del campo de velocidad en el punto. Si bien ambos espacios son tridimensionales, la aplicación de k-means produce resultados muy diferentes.

Veremos ahora cómo es que se generó esta representación alternativa, ideada durante el transcurso de esta tesis. Consideremos el sistema que generó la nube de la Figura 4.3. Al estado $S_{t}$ del sistema en un cierto instante de tiempo $t$ se lo puede describir en términos de las tres variables posicionales $\mathbf{x}_{t}=\left(x_{t}, y_{t}, z_{t}\right)$ que fueron utilizadas para producir la segmentación de la Figura 4.3. Al tratarse de un sistema dinámico, podemos considerar un espacio de características alternativo que describe a cada muestra en términos de un vector tangente a la trayectoria del sistema. Supongamos que tenemos, en adición al vector $\mathbf{x}_{t}$, un vector velocidad $\mathbf{v}_{t}$ para cada muestra. Si tomamos cada estado del sistema donde $\left\|\mathbf{v}_{t}\right\|>0$ y lo representamos en términos de un vector $\hat{\mathbf{v}}=\mathbf{v}_{t}\left\|\mathbf{v}_{t}\right\|^{-1}$ (al que llamaremos el vector flujo), obtenemos una nube de puntos alternativa $\mathcal{C}^{\prime}$ que podemos

---

[comment]: # (Page break)
segmentar con k-means para producir un resultado como el de la Figura 4.4.
![chunk1_img-20.jpeg](chunk1_img-20.jpeg)

Figura 4.4: Nube de Rössler segmentada en el espacio de vectores flujo usando K-Means. Esta segmentación no pierde el hueco central del atractor.

Esta segmentación alternativa construye segmentos alrededor del hueco, independientemente del tamaño del mismo. Intuitivamente esto se debe a que dos puntos diametralmente opuestos sobre la frontera del hueco, a pesar de ser cercanos en el espacio de vectores posición, son lejanos en el espacio de vectores flujo. En el espacio de vectores flujo una instancia del atractor de Rössler se ve como en la Figura 4.5. La intuición es que, en este espacio, dos puntos $\mathbf{x}_{a}, \mathbf{x}_{b}$ diametralmente opuestos sobre la frontera del hueco tendrán correspondientes vectores de flujo $\hat{\mathbf{v}}_{a}, \hat{\mathbf{v}}_{b}$ sobre la esfera unitaria tal que la distancia entre $\hat{\mathbf{v}}_{a}$ y $\hat{\mathbf{v}}_{b}$ es aproximadamente la máxima posible en este espacio.
![chunk1_img-21.jpeg](chunk1_img-21.jpeg)

Figura 4.5: Nube de Rössler en el espacio de vectores flujo. El hueco del atractor se vuelve relativamente grande y fácil de detectar.

Esta representación facilita la detección de huecos en los focos del atractor, alrededor de los cuales el sistema oscila. El problema con esta representación, no obstante, es que pierde la noción de distancia espacial entre muestras. En el sistema de Lorenz, por ejemplo, esto causa que conjuntos de muestras en distintas alas del atractor sean asignadas al mismo segmento, a pesar de ser espacialmente distantes (Figura 4.6).

Para aprovechar la información complementaria en cada uno de estos dos espacios, en esta tesis

---

[comment]: # (Page break)
![chunk1_img-22.jpeg](chunk1_img-22.jpeg)

Figura 4.6: Atractor de Lorenz segmentado con $k=3$, teniendo en cuenta solo el vector flujo en cada punto. Debido a la falta de información posicional en la representación, dos regiones distantes resultan asignadas al mismo segmento.
se propuso un espacio de características alternativo en $\mathbb{R}^{6}$ que se consigue construyendo el vector que tiene por primeras tres componentes las del vector posición $\mathbf{x}_{t}$, y cuenta adicionalmente con otras tres componentes correspondientes a las del vector flujo $\hat{\mathbf{v}}_{t}$. Este vector compuesto codifica por un lado la idea de distribución espacial, pero cuenta adicionalmente con información que permite distinguir más fácilmente las regiones alrededor de los focos del atractor. El atractor de Lorenz segmentado en este espacio (con previa estandarización de datos ${ }^{1}$ ) se muestra en la Figura 4.7. Lo interesante de esta segmentación es que logra "detectar" los huecos del atractor con una cantidad relativamente pequeña de segmentos.

El siguiente fragmento de código muestra cómo se construye una matriz donde cada fila es un vector de características representativo de una muestra, asumiendo que contamos con matrices x y v que representan a la nube en los espacios de vectores posición y de vectores velocidad, respectivamente:

```
# Se toman solo aquellas muestras donde la magnitud
# de la velocidad es mayor a 0.
mask = np.linalg.norm(v, axis=1) > 0.0
x, v = x[mask], v[mask]
# Se calcula el vector flujo en base al vector velocidad y a su norma.
w = np.array(list(map(lambda v_i: v_i / np.linalg.norm(v_i), v)))
# Se construye el vector compuesto concatenando posición y flujo.
return np.concatenate((x, w), axis=1)
```

En general, sea $\mathbf{s}_{t}=\left(\mathbf{x}_{t}, \mathbf{v}_{t}\right)$ con $t \in \mathbb{N} \cap[1, n]$, una sucesión de estados de un sistema dinámico muestreado a frecuencia constante $f_{s}$, donde $\mathbf{x}_{t}, \mathbf{v}_{t} \in \mathbb{R}^{m}$ y $\left\|\mathbf{v}_{t}\right\|>0 \forall t$. El vector de características

[^0]
[^0]:    ${ }^{1}$ Estandarizar un vector consiste en restar la media y dividir por el desvío estándar, componente a componente. Esta operación es un homeomorfismo de $\mathbb{R}^{m}$ sobre sí mismo y no afecta la topología del espacio; es decir, la nube original es topológicamente equivalente a la nube estandarizada.

---

[comment]: # (Page break)
se define como la version estandarizada de $\mathbf{f}_{t}=\left(\mathbf{x}_{t}, \hat{\mathbf{v}}_{t}\right)$, donde $\hat{\mathbf{v}}_{t}$ es el vector unitario en la misma dirección y sentido que $\mathbf{v}_{t}$. La segmentación de la Figura 4.7 fue producida segmentando la nube $\mathcal{X}=\left\{\mathbf{f}_{t}\right\}$ y no la nube $\mathcal{C}=\left\{\mathbf{x}_{t}\right\}$.
![chunk1_img-23.jpeg](chunk1_img-23.jpeg)

Figura 4.7: Nube de Lorenz segmentada teniendo en cuenta tanto posición como dirección del flujo.

# 4.3. El grafo de vecindad entre segmentos 

Durante el desarrollo de esta tesis se puso el foco principalmente en atractores que pueden ser aproximados por 2-complejos. Un 2-complejo puede ser construido en base a un grafo $\mathcal{G}$ (un 1-complejo) acorde al siguiente procedimiento:

1. Se selecciona mediante algún criterio un conjunto de ciclos en $\mathcal{G}$. Un ciclo es un camino cerrado definido en términos de una sucesión de vértices.
2. Por cada ciclo seleccionado se construye una 2-celda del complejo. La aristas que conforman el ciclo definen un conjunto de 1-celdas que conforman la frontera de la 2-celda. La forma concreta del interior de la 2-celda no es importante, pero se asume que se trata de una superficie acotada con frontera definida por las aristas que conforman el ciclo.

La idea fundamental detrás de los algoritmos desarrollados durante esta tesis es que, dada una segmentación adecuada de una nube de puntos $\mathcal{C}$ embebida en un espacio topológico $\mathcal{T}$, es posible construir un grafo $\mathcal{G}$ que puede ser extendido a un 2-complejo homeomórfico a $\mathcal{T}$. Esto requiere un procedimiento basado en un conjunto adecuado de criterios.

Sea $\mathcal{S}^{*}=\left\{\mathcal{S}_{i}: i=1,2, \ldots, k\right\}$, una segmentación en $k$ clusters producida por k-means aplicado a alguna nube de puntos $\mathcal{C}$ en algún espacio de características. Suponemos que cada segmento $\mathcal{S}_{i} \subset \mathcal{S}^{*}$ es una discretización de un espacio topológico $\mathcal{T}_{i}$, subespacio de aquel en el cual la nube $\mathcal{C}$ está embebida. En el límite en el que la densidad de muestras tiende a infinito, y asumiéndola uniforme, podemos intuir que $\mathcal{S}_{i}$ se asemeja cada vez más a $\mathcal{T}_{i}$ (Figura 4.8). En otras palabras, podemos decir que $\mathcal{S}^{*}$ aproxima una partición del espacio topológico $\mathcal{T}^{*}=\bigcup_{i=1}^{k} \mathcal{T}_{i}$ subyacente a la nube de puntos. El algoritmo ideal sería entonces uno que por cada segmento $\mathcal{S}_{i}$ construye una 2-celda cuya frontera se asemeja a la de $\mathcal{S}_{i}$. Durante el trabajo se logró producir una variedad de

---

[comment]: # (Page break)
![chunk1_img-24.jpeg](chunk1_img-24.jpeg)

Figura 4.8: Asumimos que, a medida que la densidad de muestras aumenta, $\mathcal{S}_{i}$ tiende al espacio topológico subyacente $\mathcal{T}_{i}$.
![chunk1_img-25.jpeg](chunk1_img-25.jpeg)

Figura 4.9: Un grafo construido segmentando una nube de puntos y conectando ciertos "puntos de interés" a través de aristas embebidas en ciertas "regiones de interés". No se logró que esto funcione en la juntura. El grafo se dibujó superpuesto sobre la nube de puntos.
ideas que permiten producir 1-esqueletos "a simple vista" satisfactorios, como el que se muestra en la Figura 4.9. El algoritmo que generó el grafo de la Figura 4.9 produce resultados satisfactorios en regiones homeomórficas al plano, pero se encontró significativamente más difícil desarrollar un método que produzca resultados consistentes en las junturas, las regiones donde las ramas confluyen. Los resultados más consistentes se obtuvieron con una estrategia alternativa: por cada segmento se produce una 0 -celda en vez de una 2 -celda, y las 0 -celdas se conectan acorde a un criterio de vecindad $\nu: \mathcal{S}^{*} \times \mathcal{S}^{*} \rightarrow\{0,1\}$ tal que $\nu\left(\mathcal{S}_{i}, \mathcal{S}_{j}\right)=1$ si los segmentos $\mathcal{S}_{i}$ y $\mathcal{S}_{j}$ son vecinos (caso en el que sus correspondientes 0 -celdas se conectan mediante una 1-celda), o 0 en caso contrario. Con un criterio $\nu$ adecuado y una cantidad suficientemente grande de segmentos es posible construir grafos como el que se muestra en la Figura 4.10. Este grafo es lo que llamaremos un grafo de vecindad entre segmentos. Se encontró empíricamente que, dada una segmentación adecuada, el 1-complejo definido por el grafo de vecindad entre segmentos puede ser extendido a un 2-complejo en base al cual la etapa final de BraMAH calcula los invariantes topológicos correctos, correspondientes al atractor que la nube aproxima. Los detalles del criterio de vecindad adoptado concretamente se detallan en la Sección 5.3.

---

[comment]: # (Page break)
![chunk1_img-26.jpeg](chunk1_img-26.jpeg)

Figura 4.10: Un grafo de vecindad entre segmentos construido en base a una nube producida tomando muestras de la trayectoria de un sistema de Lorenz. Seleccionando los ciclos adecuados para conformar 2-celdas se puede producir un complejo homeomórfico al atractor de Lorenz. El grafo se dibujó superpuesto sobre la nube de puntos.

---

[comment]: # (Page break)
# Sección 5 

## Construcción del complejo celular

A continuación se describe en detalle el algoritmo utilizado para producir un complejo celular en base a una nube de puntos. La Figura 5.1 describe el flujo de datos a través de las distintas etapas de ejecución del algoritmo. Este proceso es el que se utilizó para construir el complejo de la Figura 5.7 en base a una nube obtenida muestreando una instancia del sistema de Rössler. Las etapas indicadas en el diagrama son las siguientes:

- Preparación del conjunto de datos (Secciones 5.1 y 5.2). Este paso debe ser realizado por el usuario y requiere producir un conjunto de datos con las características indicadas.
- Segmentación (Sección 5.3). En este paso se utilizan las técnicas detalladas en la Sección 4 para segmentar la nube de puntos provista por el usuario.
- Construcción de un 1-complejo primitivo (Sección 5.4). En este paso se construye un grafo que codifica la estructura de un 1-complejo celular. Este grafo será luego extendido a un 2-complejo identificando las 2-celdas o "caras" del complejo.
- Identificación de huecos topológicos (Secciones 5.5 y 5.6). En este paso se identifican aquellos ciclos del 1-complejo producido en el paso anterior que representan 2-celdas y aquellos que son huecos topológicos. Esta es la información necesaria para extender el 1-complejo a un 2-complejo.
- Orientación de celdas (Sección 5.7). En este paso se le asigna a cada 2-celda una de dos posibles orientaciones para definir un complejo orientado.

Se implementó también un procedimiento de optimización que busca mejorar los resultados producidos por el algoritmo inicial. Este procedimiento se describe en la Sección 5.8.

### 5.1. Generación del conjunto de datos

El código concreto utilizado para producir el conjunto de datos de entrada se muestra a continuación:

---

[comment]: # (Page break)
![chunk1_img-27.jpeg](chunk1_img-27.jpeg)

Figura 5.1: Construcción de un complejo celular en base a una nube de puntos.

```
import numpy as np
from scipy.integrate import solve_ivp
def rossler(t, r, a=0.43295, b=2.08, c=4.0):
    x, y, z = r
    dzdt = -y - z
    dydt = x + a * y
    dzdt = b + z + (x - c)
    return [dzdt, dydt, dzdt]
# Se definen las propiedades de la integración.
t_span, X0 = np.linspace(0, 1000, 20000), [0.0, -2.0, 0.5]
# Se integra numéricamente.
solved = solve_ivp(rossler, [t_span[0], t_span[-1]], X0, t_eval=t_span)
# Se construye la matriz donde cada fila es un vector posición.
x = np.vstack([solved.y[0], solved.y[1], solved.y[2]]).T
```

El código hace uso de los módulos numpy y scipy (por solve_ivp). Este código produce una secuencia de estados discretos ordenados cronológicamente que en la Figura 5.2 se muestra como una curva continua.

---

[comment]: # (Page break)
![chunk1_img-28.jpeg](chunk1_img-28.jpeg)

Figura 5.2: Trayectoria que describe el atractor de Rössler, producida integrando numéricamente el sistema de Rössler con parámetros $(a, b, c)=(0.43295,2.08,4.0)$ y condiciones iniciales $\left(x_{0}, y_{0}, z_{0}\right)=(0,-2,0.5)$, desde $t=0$ hasta $t=1,000$ en 19,999 pasos, usando la implementación RK4S del método de Runge-Kutta en la biblioteca SciPy 1.10.0.

# 5.2. Supuestos sobre los datos de entrada 

Sea $\mathcal{C}=\left\{\mathbf{x}_{t}\right\} \subset \mathbb{R}^{m}$, con $t \in 1,2, \ldots, n$ una nube de puntos y $\mathbf{x}_{t}$ una sucesión tal que para todo $t<n, \mathbf{x}_{t}$ precede inmediatamente a $\mathbf{x}_{t+1}$ en orden cronológico. Es decir, suponemos que la sucesión $\mathbf{x}_{t}$ describe la discretización de una trayectoria en el espacio de estados de un sistema dinámico $S$. Suponemos adicionalmente que la nube $\mathcal{C}$ aproxima al conjunto de puntos de un espacio topológico subyacente $\mathcal{T}$, y que existe un 2-complejo homeomórfico a $\mathcal{T}$. Es decir, suponemos que $\mathcal{T}$ puede ser representado adecuadamente por un 2-complejo. El objetivo del algoritmo será producir un 2-complejo orientado $\mathcal{K}$, homeomórfico a $\mathcal{T}$.

### 5.3. Segmentación

El primer paso del algoritmo consiste en segmentar la nube como se describió en la Sección 4.2. Si no se cuenta adicionalmente con un conjunto de vectores velocidad $\mathbf{v}_{t}$ para cada muestra, es posible estimar $\mathbf{v}_{t}$ tomando la primera diferencia de la sucesión $\mathbf{x}_{t}$ tal que

$$
\mathbf{v}_{t}=\mathbf{x}_{t+1}-\mathbf{x}_{t}, \quad t<n
$$

Esto requiere descartar la última muestra $\mathbf{x}_{n}$, o definir algún criterio de aproximación que permita definir $\mathbf{v}_{n}$. Un posible criterio consiste en suponer $\mathbf{v}_{n} \approx \mathbf{v}_{n-1}$, asumiendo la suavidad del campo de velocidad. El siguiente fragmento de código muestra esta última aproximación:

```
# Se toman primeras diferencias del vector posición.
v = np.diff(x, axis=0)
# Se aproxima v[n] como v[n-1].
v = np.append(v, [v[-1]], axis=0)
```

---

[comment]: # (Page break)
Procedemos luego definiendo una nueva nube $\mathcal{X}$ de $n^{\prime}$ elementos en un espacio $2 m$ dimensional. Sea $\mathbf{s}_{t}^{\prime}=\left(\mathbf{x}_{t}^{\prime}, \mathbf{v}_{t}^{\prime}\right)$ una subsucesión de $\mathbf{s}_{t}=\left(\mathbf{x}_{t}, \mathbf{v}_{t}\right)$ tal que $\mathbf{s}_{t}^{\prime}$ contiene los elementos de $\mathbf{s}_{t}$ donde $\left\|\mathbf{v}_{t}\right\|>0$. Existe un elemento en $\mathcal{X}$ por cada muestra $\mathbf{s}_{t}^{\prime}$. Dado un vector $\mathbf{f}_{i} \in \mathcal{X}$, sus componentes son tal como se describe en la Sección 4.2: se concatenan los vectores $\mathbf{x}_{i}^{\prime}$ y $\dot{\mathbf{v}}_{i}^{\prime}$, y el resultado se estandariza restando la media y dividiendo por el desvío estándar componente a componente. El siguiente fragmento de código describe el procedimiento completo:

```
# Se toman primeras diferencias del vector posición.
v = np.diff(x, axis=0)
# Se aproxima v[n] como v[n-1].
v = np.append(v, [v[-1]], axis=0)
# Se toman solo aquellas muestras donde la magnitud
# de la velocidad es mayor a 0.
mask = np.linalg.norm(v, axis=1) > 0.0
x, v = x[mask], v[mask]
# Se calcula el vector flujo en base al vector velocidad y a su norma.
w = np.array(list(map(lambda v_i: v_i / np.linalg.norm(v_i), v)))
# Se construye el vector compuesto concatenando posición y flujo.
f = np.concatenate((x, w), axis=1)
# Se aplica el procedimiento de estandarización usando
# el módulo preprocessing de scikit learn.
scaler = preprocessing. StandardScaler().fit(f)
scaled = scaler.transform(f)
# Se segmenta el conjunto de datos usando el módulo
# clustering the scikit learn.
kmeans = clustering.KMeans(init="random", n_clusters=k, max_iter=max_iter)
kmeans.fit(scaled)
# Se obtiene la asignación de segmentos a muestras.
coloring = kmeans.predict(scaled)
```

A la nube $\mathcal{X}$ resultante se la segmenta usando k-means con un valor de $k$-definido por el usuario. Esto produce un conjunto $\mathcal{S}^{*}=\left\{\mathcal{S}_{j}\right\}$ de $q \leq k$ segmentos. En base a esta segmentación se le asigna a cada muestra $\mathbf{s}_{t}^{\prime}$ un color $c_{i}$. Decimos que $\mathbf{s}_{t}^{\prime}$ es de color $c_{i}$ si $\mathbf{f}_{i} \in \mathcal{S}_{c_{i}}$. La segmentación produce segmentos disjuntos, por lo que cada muestra es de un único color. Llamamos $c_{i}$ al color de la muestra $\mathbf{s}_{t}^{\prime}$, y llamamos $\mathcal{L}_{c}$ al conjunto de todas las muestras de color $c$. A cada $\mathcal{L}_{c}$ le corresponde un segmento $\mathcal{S}_{c}$, y al conjunto $\bigcup_{c} \mathcal{L}_{c}$ de todos $\operatorname{los} \mathcal{L}_{c}$ lo llamamos $\mathcal{L}^{*}$ en analogía con el conjunto $\mathcal{S}^{*}$. Los conjuntos $\mathcal{L}_{c}$ y $\mathcal{S}_{c}$ son equivalentes en el sentido de que refieren al mismo conjunto de entidades, siendo $\mathcal{L}_{c}$ un conjunto de muestras $\left\{\mathbf{s}_{t}^{\prime}: c_{i}=c\right\}$ y $\mathcal{S}_{c}$ un conjunto de correspondientes vectores de características $\mathbf{f} \in \mathcal{X}$. Es decir, podemos definir una asignación $\mathcal{L}_{c} \mapsto \mathcal{S}_{c}$ uno a uno e inversible. Para referirnos a un segmento o cluster particular, se usarán $\mathcal{L}_{c}$ o $\mathcal{S}_{c}$ en forma intercambiable. El producto de segmentar el atractor de Rössler de la figura 5.2 aplicando esta técnica se muestra en la figura 5.3.

---

[comment]: # (Page break)
![chunk1_img-29.jpeg](chunk1_img-29.jpeg)

Figura 5.3: Atractor de Rössler segmentado usando KMeans con $k=40$ y un máximo de 100 iteraciones. La implementación usada es la que provee el paquete scikit-learn 1.2.1, especificando una semilla random_seed=1. Todos los restantes parámetros son aquellos que la implementación establece por defecto.

# 5.4. Construcción del 1-complejo 

El siguiente paso consiste en construir un 1-complejo que pueda ser luego extendido al 2-complejo buscado. Se procede como se describió en la Sección 4.3. Se genera una 0 -celda por segmento del complejo, y las 0 -celdas se conectan acorde a un criterio de vecindad entre segmentos.

El criterio adoptado para determinar si dos segmentos son vecinos es el siguiente. Decimos que dos segmentos $\mathcal{L}_{a}, \mathcal{L}_{b}$ son vecinos si existe flujo desde $\mathcal{L}_{a}$ hacia $\mathcal{L}_{b}$, o desde $\mathcal{L}_{b}$ hacia $\mathcal{L}_{a}$. Sea $\Phi: \mathcal{L}^{*} \times \mathcal{L}^{*} \rightarrow\{0,1\}$ una función booleana tal que $\Phi\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)$ es 1 si existe flujo desde $\mathcal{L}_{a}$ hacia $\mathcal{L}_{b}$, o 0 en caso contrario. Para definir $\Phi$ usamos el ordenamiento cronológico de $\mathbf{s}_{i}^{\prime}$. Definimos primero un criterio inicial $\phi: \mathcal{L}^{*} \times \mathcal{L}^{*} \rightarrow\{0,1\}$ tal que

$$
\phi\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)=1 \quad \Leftrightarrow \quad \exists \mathbf{s}_{i}^{\prime} \in \mathcal{L}_{a} \text { tal que } \mathbf{s}_{i+1}^{\prime} \in \mathcal{L}_{b}
$$

Es decir, $\phi\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)$ es 1 si y solo si existe alguna muestra $\mathbf{s}_{i}^{\prime}$ de color $a$ tal que la muestra $\mathbf{s}_{i+1}^{\prime}$, inmediatamente posterior a $\mathbf{s}_{i}^{\prime}$ en orden cronológico, es de color $b$.

La función $\phi$ define un criterio inicial por el cual decimos que hay flujo desde $\mathcal{L}_{a}$ hacia $\mathcal{L}_{b}$ si existe alguna línea de flujo que ingresa a $\mathcal{L}_{b}$ desde $\mathcal{L}_{a}$. Para la implementación se adoptó un criterio que construye sobre $\phi$ para proveer mayor robustez ante posibles defectos de la segmentación. Sea $p\left(\mathcal{L}_{a}\right)$ la partición de $\mathcal{L}_{a}$ compuesta por todos los subconjuntos de tamaño 1 de $\mathcal{L}_{a}$. Definimos el criterio $\Phi_{\tau}: \mathcal{L}^{*} \times \mathcal{L}^{*} \rightarrow\{0,1\}$ tal que

$$
\Phi_{\tau}\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)=1 \quad \Leftrightarrow \quad \sum_{s \in p\left(\mathcal{L}_{a}\right)} \phi\left(s, \mathcal{L}_{b}\right) \geq \tau
$$

donde $\tau$ es un parámetro que define un umbral. Es decir, existe flujo desde $\mathcal{L}_{a}$ hacia $\mathcal{L}_{b}$ si la cantidad de líneas de flujo que ingresan a $\mathcal{L}_{b}$ desde $\mathcal{L}_{a}$ es mayor o igual a un cierto umbral $\tau$. Este es el criterio adoptado para la implementación, siendo el parámetro $\tau$ configurable por el usuario. Un $\tau$ mayor implica que el algoritmo requiere más evidencia para aceptar la existencia de flujo entre dos

---

[comment]: # (Page break)
segmentos. En la implementación el valor por defecto de $\tau$ es 2 . Se encontró experimentalmente que este valor produce los mejores resultados para los casos de prueba evaluados, mitigando conexiones espurias sin impactar en la topología resultante.

El siguiente fragmento de código detalla el procedimiento:

```
# Se instancia un mapa que asigna a cada color el conjunto
# de sus colores vecinos.
flow_map = defaultdict(set)
# Se itera a través de los respectivos colores de
# todas las muestras.
for i in range(len(coloring) - n):
    # Se obtiene el color de la i-ésima muestra.
    c_i = coloring[i]
    # Se asegura que existe al menos un conjunto vacío para el color de la
    # muestra actual c_i.
    flow_map[c_i] = set() if c_i not in flow_map else flow_map[c_i]
    # Se obtienen los respectivos colores de las siguientes n muestras
    # de la nube.
    upcoming_colors = set([coloring[i + j] for j in range(1, n + 1)])
    # Se verifica que el conjunto de los colores de las siguientes n muestras
    # contenga un único color. Si es este el caso, las siguientes n muestras
    # pertenecen todas al mismo segmento.
    if len(upcoming_colors) > 1:
        continue
    # Se verifica que el color de las siguientes n muestras
    # sea efectivamente distinto al de la muestra actual.
    c_j = upcoming_colors.pop()
    if upcoming_color_j == c_i:
        continue
    # Se cumplen las condiciones para indicar que existe flujo
    # desde el cluster con color c_i a aquel con color c_j.
    flow_map[c_i].add(c_j)
```

Tras la ejecución del algoritmo, el mapa flow_map contendrá una descripción del flujo entre segmentos tal que flow_map[c] es el conjunto de los colores de segmentos hacia los cuales existe flujo desde el segmento de color c.

Definida la noción de vecindad entre segmentos, definimos el 1-complejo $\mathcal{K}^{1}$ cuyas 0 -celdas y 1-celdas son respectivamente los vértices y aristas de un grafo $\mathcal{G}=(\mathcal{V}, \mathcal{E})$ que tiene un vértice $v_{c}$ por cada segmento $\mathcal{L}_{c}$, y dos vértices $v_{a}, v_{b}$ se conectan mediante una arista si $\Phi_{\tau}\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)=1$ o si $\Phi_{\tau}\left(\mathcal{L}_{b}, \mathcal{L}_{a}\right)=1$ (Figura 5.4).

Sea $\mathcal{C}^{\prime}$ el conjunto de todos $\operatorname{los} \mathbf{x}_{i}^{\prime}$ de la sucesión $\mathbf{s}_{i}^{\prime}=\left(\mathbf{x}_{i}^{\prime}, \mathbf{v}_{i}^{\prime}\right)$. Para dotar al grafo $\mathcal{G}$ de información geométrica, le asignamos a cada vértice $v_{c} \in \mathcal{V}$ una posición $\mathbf{r}\left(v_{c}\right)$ tal que $\mathbf{r}\left(v_{c}\right)$ es el punto $\mathbf{x} \in \mathcal{C}^{\prime}$ más cercano al centroide del conjunto $\left\{\mathbf{x}_{i}^{\prime} \in \mathcal{C}^{\prime}: \mathbf{s}_{i}^{\prime} \in \mathcal{L}_{c}\right\}$. Es decir, de cada cluster $\mathcal{L}_{c}$ se calcula su centroide en el espacio de vectores posición, y el punto de la nube $\mathcal{C}^{\prime}$ más cercano a dicho centroide será la posición representativa del vértice $v_{c}$.

---

[comment]: # (Page break)
![chunk1_img-30.jpeg](chunk1_img-30.jpeg)

Figura 5.4: Grafo construido asignando a cada segmento un punto representativo, y uniendo puntos mediante aristas acorde a un criterio de vecindad entre segmentos. Los vértices y aristas del grafo definen respectivamente las 0-celdas y 1-celdas del complejo $\mathcal{K}^{1}$.

# 5.5. Identificación de las 2-fronteras candidatas 

Habiendo construido el grafo $\mathcal{G}$, el siguiente paso consiste en identificar los ciclos que serán las fronteras de las 2-celdas del complejo $\mathcal{K}$. Lo que buscamos es un concepto análogo al de "cara" en grafos planos. El problema es que en grafos no planos, tal concepto no es tan simple de definir. Como heurística se adoptó el concepto de la base mínima de ciclos [27]. Una base de ciclos del grafo $\mathcal{G}$ es un conjunto de ciclos en $\mathcal{G}$ que, bajo la operación de diferencia simétrica (la unión menos la intersección) de aristas, permite construir cualquier ciclo del grafo (Figura 5.5). Si a cada ciclo le asignamos un peso (determinado acorde a alguna métrica arbitraria), la base de ciclos es mínima respecto a la métrica si la sumatoria de los pesos de los ciclos en la base es mínima en el espacio de todas las posibles bases de ciclos de $\mathcal{G}$. Se encontró experimentalmente que una base mínima de ciclos de $\mathcal{G}$, tomando como peso de un ciclo a la sumatoria de la longitud geométrica de las aristas que componen al ciclo, es un buen candidato para el conjunto de las caras del grafo $\mathcal{G}$.

Para definir el peso de cada ciclo usamos la información geométrica en los vértices. Dada una arista $e \in \mathcal{E}$ que une dos vértices $v_{a}, v_{b} \in \mathcal{V}$, le asignamos a $e$ una longitud geométrica $\ell(e)$ tal que

$$
\ell(e)=\left\|\mathbf{r}\left(v_{a}\right)-\mathbf{r}\left(v_{b}\right)\right\|
$$

Luego, dado un ciclo $C$ definido como un conjunto ordenado de aristas, definimos el peso $w(C)$ del ciclo $C$ como la suma

$$
w(C)=\sum_{e \in C} \ell(e)
$$

La base mínima de ciclos de $\mathcal{G}$ se calcula minimizando el peso total de la base. Esto produce un conjunto de ciclos $\mathcal{B}$ que son además potenciales 2-fronteras del complejo celular $\mathcal{K}$.

En esta instancia se imponen ciertas restricciones para mitigar posibles defectos introducidos por una segmentación inicial inadecuada. Sea $C \in \mathcal{B}$ la frontera de una 2-celda del complejo $\mathcal{K}$,

---

[comment]: # (Page break)
![chunk1_img-31.jpeg](chunk1_img-31.jpeg)

Figura 5.5: Ejemplo de un grafo generado por dos ciclos A $y$ B. Al ciclo A lo define la secuencia de vértices 1,2,5,6,1. Al ciclo B lo define la secuencia 2,3,4,5,2. El grafo contiene adicionalmente un ciclo 1,2,3,4,5,6,1 generado por la unión entre los conjuntos de aristas de A $y$ B, quitando la intersección entre los conjuntos de aristas de A $y$ B.
compuesta por los vértices $\mathcal{V}_{C} \subseteq \mathcal{V}$ interconectados por las aristas $\mathcal{E}_{C} \subseteq \mathcal{E}$. Es decir, $\mathcal{G}_{C}=\left(\mathcal{V}_{C}, \mathcal{E}_{C}\right)$ es un subgrafo de $\mathcal{G}$ compuesto por los vértices y las aristas que conforman al ciclo $C$. Imponemos la restricción de que todos los vértices $\mathcal{G}_{C}$ deben ser de grado dos; es decir, unidas a cada vértice debe haber dos y solo dos aristas. De no cumplirse esta condición para todo $C \in \mathcal{B}$, se reintenta el procedimiento con una nueva segmentación y una nueva base de ciclos hasta construir un grafo que la cumpla. Esto último implica lo siguiente:

- No hay garantía de que el algoritmo termine para cualquier conjunto de parámetros. Esta complicación no se observó para ninguno de los casos de prueba evaluados.
- Se hacen supuestos sobre la dimensionalidad de las celdas. En otras palabras, no se trata de un algoritmo general capaz de tratar cualquier tipo de espacio.

Estas son quizás las limitaciones más evidentes del algoritmo. En las pruebas se encontró, no obstante, que la ejecución termina consistentemente y produce resultados satisfactorios, para los casos evaluados.

Como paso adicional, para cada ciclo los vértices se ordenan en una lista tal que el orden de los vértices defina un recorrido alrededor del ciclo. Es decir, si el vértice $v_{b}$ es inmediatamente posterior en orden a un vértice $v_{a}$, un recorrido alrededor del ciclo implica en algún punto del trayecto recorrer una arista que une los vértices $v_{a}$ y $v_{b}$.

# 5.6. Identificación de las 2-celdas 

Contamos ahora con un conjunto $\mathcal{B}$ de ciclos candidatos a ser las fronteras de las 2-celdas de $\mathcal{K}$. El siguiente paso consiste en distinguir cuáles de estas fronteras encierran en realidad huecos topológicos del espacio subyacente. Las 2-celdas del complejo quedarán automáticamente definidas como el interior de las 2-fronteras que no encierren huecos.

Sea $\mathcal{F}=\left(\mathcal{V}, \mathcal{E}_{\mathcal{F}}\right)$ el grafo compuesto por los vértices de $\mathcal{G}$, conectados por aristas dirigidas acorde al siguiente criterio: dos vértices $v_{a}, v_{b} \in \mathcal{V}$ se conectan mediante una arista dirigida $\left(v_{a}, v_{b}\right) \in \mathcal{E}_{\mathcal{F}}$ si $\varphi_{\tau}\left(\mathcal{L}_{a}, \mathcal{L}_{b}\right)=1$. Es decir, $\mathcal{F}$ es una versión dirigida del grafo $\mathcal{G}$, donde dos vértices $v_{a}, v_{b}$ se conectan

---

[comment]: # (Page break)
mediante una arista dirigida $\left(v_{a}, v_{b}\right)$ si existe flujo desde el segmento $\mathcal{L}_{a}$ hacia el segmento $\mathcal{L}_{b}$. En otras palabras, $\mathcal{F}$ es un grafo de flujo entre segmentos (Figura 5.6).
![chunk1_img-32.jpeg](chunk1_img-32.jpeg)

Figura 5.6: Grafo de flujo construido en base a la segmentación del atractor de Rössler de la figura 5.3. El flujo alrededor del foco del atractor define un ciclo orientado cerrado.

Dado un ciclo $C \in \mathcal{B}$, adoptamos como criterio para determinar si $C$ encierra un hueco el que se describe a continuación: Sea $v^{1}, v^{2}, \ldots, v^{\kappa}$ una sucesión de $\kappa$ vértices $v^{i} \in \mathcal{V}$ que describe un recorrido simple del ciclo $C$. Decimos que $C$ es un hueco si la secuencia de vértices $v^{i}$ describe un ciclo en el grafo $\mathcal{F}$, teniendo en cuenta la dirección de las aristas de $\mathcal{F}$. Todos aquellos ciclos que no encierran huecos conforman un conjunto $\mathcal{B}^{\prime} \subset \mathcal{B}$.

# 5.7. Orientación de las celdas 

El paso final en la construcción del complejo consiste en orientar sus 2-celdas en forma consistente. La orientación de una celda está implícita en el orden en el que se especifican los vértices que conforman su frontera. Dada una 2-celda $A$, la especificamos en términos de una sucesión de vértices $v_{i}^{A}$ que describe un ciclo simple $C \in \mathcal{B}^{\prime}$ sobre su frontera, en algún sentido. Dada cualquier 2-celda $B$ distinta de $A$ tal que $A$ y $B$ comparten al menos un par de vértices, todo par de vértices $v_{a}, v_{b}$ compartido por las celdas $A$ y $B$ debe cumplir la siguiente implicación:

$$
\exists i \quad v_{i}^{A}=v_{a}, \quad v_{i+1}^{A}=v_{b} \quad \Rightarrow \quad \exists j \quad v_{j}^{B}=v_{b}, \quad v_{j+1}^{B}=v_{a}
$$

Es decir, cualquier arista compartida por dos celdas fronterizas que es recorrida en un sentido en una celda debe ser recorrida en sentido opuesto en la otra. Una celda queda definida entonces por un ciclo orientado $C \in \mathcal{B}^{\prime}$, y decimos que dos ciclos son vecinos si comparten al menos dos vértices.

El algoritmo que orienta las celdas del complejo procede de la siguiente forma:

1. El algoritmo termina si todas las celdas han sido orientadas. Mientras haya alguna celda sin orientar, se procede con los siguientes pasos.
2. Se toma arbitrariamente un ciclo $C \in \mathcal{B}^{\prime}$ que no haya sido todavía orientado y se le asigna una orientación arbitraria. El ciclo con su orientación define una 2-celda del complejo orientado. Se marca el ciclo $C$ como orientado.

---

[comment]: # (Page break)
3. Se propaga la orientación de $C$ hacia sus ciclos vecinos como se describe a continuación. Se identifican los ciclos $\left\{X_{i}\right\} \subset \mathcal{B}^{\prime}$ todavía no orientados que comparten aristas con $C$, y se define el orden de los vértices de cada $X_{i}$ tal de mantener la consistencia en la orientación con $C$. A cada ciclo $X_{i}$ se lo marca como orientado luego de definir su orientación.
4. Se repite el paso 3 en forma recursiva por cada celda $X_{i}$ con $C \leftarrow X_{i}$. Si no hay $X_{i}$ con algún ciclo vecino todavía no orientado, se retorna al paso 1 .

El resultado de este paso es entonces un complejo orientado.

# 5.8. Optimización aleatoria 

El complejo resultante tras aplicar la secuencia de pasos descrita hasta ahora se muestra en la figura 5.7. La segmentación detrás de este complejo fue generada definiendo una semilla con valor 0 para la implementación de k-means de Scikit-Learn, produciendo un resultado determinístico. En general, no obstante, no definir una semilla en forma explícita causará que el conjunto inicial de centroides sea seleccionado en forma aleatoria; es decir, el complejo resultante puede variar de ejecución a ejecución.
![chunk1_img-33.jpeg](chunk1_img-33.jpeg)

Figura 5.7: Complejo celular construido en base a la nube de puntos de la figura 5.3. Las 2-celdas se muestran en color azul. Los puntos grises indican los centros de las celdas. La 1-celda indicada en rojo es compartida por tres o más 2-celdas y se encuentra sobre una línea de juntura.

Para los conjuntos de datos evaluados, en hardware comercial accesible ${ }^{1}$, la ejecución del algo-

[^0]
[^0]:    ${ }^{1}$ MacBook Air genérica con procesador Apple M1, 8GB de RAM, MacOS 13.

---

[comment]: # (Page break)
ritmo termina en pocos segundos. Para mejorar la calidad del resultado se favoreció entonces una extensión estilo Monte Carlo que ejecuta el proceso repetidas veces y a cada complejo resultante $\mathcal{K}_{i}$ lo califica acorde a una función de aptitud $\alpha: \mathcal{K}^{*} \rightarrow \mathbb{R}$, donde $\mathcal{K}^{*}$ es el conjunto de todos los complejos celulares que puede producir el algoritmo. Dada una cantidad predefinida $\iota$ de iteraciones que producen cada una un complejo $\mathcal{K}_{i}$, con $i \in \mathbb{N} \cap[1, \iota]$, el producto final del proceso de optimización es un complejo

$$
\mathcal{K}_{x}, \quad x=\arg \max _{i} \alpha\left(\mathcal{K}_{i}\right), \quad i \in \mathbb{N} \cap[1, \iota]
$$

Como criterio de aptitud se utilizó $-J_{i}$, donde $J_{i}$ es la cantidad de aristas utilizadas para representar las denominadas líneas de juntura del complejo $\mathcal{K}_{i}$. Sea $\mathcal{K}$ un 2-complejo definido en términos de un conjunto de ciclos orientados $\mathcal{B}^{\prime}$ de un grafo $\mathcal{G}=(\mathcal{V}, \mathcal{E})$. Una arista $e \in \mathcal{E}$ sobre una línea de juntura es una arista compartida por 3 o más ciclos de $\mathcal{B}^{\prime}$. Es decir, $e$ es una arista sobre la cual el complejo no es una superficie. Sea $\jmath: \mathcal{E} \rightarrow\{0,1\}$ una función booleana que es 1 si $e$ es una arista compartida por tres o más ciclos en $\mathcal{B}^{\prime}$, o 0 en caso contrario. Definimos la aptitud del complejo $\mathcal{K}$ como

$$
\alpha(\mathcal{K})=-\sum_{e \in \mathcal{E}} \jmath(e)
$$

Es decir, se favorecen (tienen aptitud más alta) los complejos basados en grafos que requieren menos aristas para representar las líneas de juntura. Independientemente de que tales complejos pueden resultar más parsimoniosos, este criterio es en gran parte arbitrario y la implementación detallada en la Sección 7 es parametrizable con un criterio a elección del usuario.

Nótese que el procedimiento detallado no necesariamente identificará un óptimo global en $\mathcal{K}^{*}$, sino el mejor entre los $\iota$ complejos generados, acorde al criterio de aptitud seleccionado. Es decir, dada una muestra aleatoria de $\iota$ complejos de $\mathcal{K}^{*}$, seleccionamos aquel que consideramos más adecuado acorde al criterio de calidad $\alpha(\mathcal{K})$. Si bien definir el complejo buscado como un óptimo global de alguna función adecuada es una posible forma de abordar la construcción de complejos celulares en base a datos, la complejidad de la tarea excede el alcance de esta tesis.

---

[comment]: # (Page break)
# Sección 6 

## Construcción del templex

### 6.1. Extensión de un complejo a un templex

En la Figura 6.1 se muestra una representación gráfica del digrafo (grafo dirigido) de un templex construido en base al complejo celular de la Figura 5.7 y a la nube de puntos utilizada para generarlo. La representación incluye por un lado el complejo celular con sus correspondientes líneas de juntura y sobre él, superpuesto, se encuentra el digrafo de flujo que describe las transiciones válidas entre 2-celdas a través de 1-celdas en el espacio de fases. El digrafo fue generado mediante el algoritmo que se describe a continuación.
![chunk1_img-34.jpeg](chunk1_img-34.jpeg)

Figura 6.1: Representación gráfica de un templex construido en base a un atractor de Rössler. El gráfico consiste en el complejo celular de la Figura 5.7 sobre el cual se superpone un grafo dirigido o digrafo que describe la dirección del flujo a través de las 1-celdas del complejo.

---

[comment]: # (Page break)
Sea $\mathcal{K}$ un complejo celular construido en base a una nube de puntos $\mathcal{C}=\left\{\mathbf{x}_{t}\right\} \subset \mathbb{R}^{m}, t=$ $1,2, \ldots, n$, donde $\mathbf{x}_{t}$ discretiza una trayectoria en el espacio de fases. Cada celda del complejo celular $\mathcal{K}$ tiene asociada información geométrica: una 0 -celda $P$ tiene asociada una posición $\mathbf{r}(P) \in \mathcal{C}$, y cada 2-celda $\sigma$ tiene asociado un punto representativo $\mathbf{r}(\sigma) \in \mathcal{C}$.

Sea $j: \mathcal{K}_{1} \rightarrow \mathcal{P}\left(\mathcal{K}_{2}\right)$, con $\mathcal{K}_{1}$ el conjunto de 1-celdas y $\mathcal{P}\left(\mathcal{K}_{2}\right)$ el conjunto de todos los subconjuntos de 2-celdas de $\mathcal{K}$. Definimos como $j(e)$ al conjunto de todas las 2-celdas de $\mathcal{K}$ que tienen a la 1-celda $e \in \mathcal{K}_{1}$ en su frontera. El algoritmo requiere definir un criterio para asignar a cada posible $j(e)$ una 2-celda $\sigma^{*} \in j(e)$, a la que llamaremos receptora de flujo. Esta celda es la que recibe el flujo proveniente de las restantes 2-celdas en $j(e)$. Este criterio deberá usar en general información contenida en la nube $\mathcal{C}$, así como información geométrica asociada al complejo $\mathcal{K}$. El criterio concreto adoptado para la implementación se describirá en la Sección 6.2.

El algoritmo de extensión produce un grafo dirigido $\mathcal{F}_{\Phi}=\left\langle\mathcal{V}_{\Phi}, \mathcal{E}_{\Phi}\right\rangle$ que junto con $\mathcal{K}$ define un simplex. En el conjunto $\mathcal{V}_{\Phi}$ existirá un vértice $v_{\sigma}$ por cada 2-celda $\sigma$ en $\mathcal{K}_{2}$, y cada vértice $v_{\sigma}$ tiene asociada una posición $\mathbf{r}\left(v_{\sigma}\right)=\mathbf{r}(\sigma)$. El algoritmo comienza ejecutando la siguiente secuencia de pasos por cada 1-celda $e$ del complejo $\mathcal{K}$ :

1. Se determina el conjunto $j(e)$.
2. Se determina la celda $\sigma^{*} \in j(e)$ receptora de flujo en $j(e)$, acorde a un criterio a definir. Este criterio hará uso de la información geométrica, topológica y temporal en $\mathcal{K}$ y $\mathcal{C}$.
3. Por cada $\sigma \in j(e)$, con $\sigma \neq \sigma^{*}$, se agrega a $\mathcal{E}_{\Phi}$ una arista dirigida $\left(v_{\sigma}, v_{\sigma^{*}}\right)$.

La dificultad del algoritmo está en la definición de un criterio adecuado para calcular el flujo entre celdas. Los detalles del criterio adoptado se describen a continuación en la Sección 6.2.

# 6.2. Criterio de orientación de flujo 

Dada la 1-celda $e \in \mathcal{K}_{1}$ y el conjunto $j(e)$ de 2-celdas que tienen a $e$ en su frontera, debemos definir un criterio que nos permita definir la dirección del flujo a través de $e$. Suponemos en principio que la granularidad de las celdas es tal que la dirección del vector flujo es aproximadamente uniforme sobre $e$. El criterio asume que es correcto afirmar la existencia de una única celda $\sigma^{*} \in j(e)$ hacia la cual confluye el flujo proveniente de las restantes celdas en $j(e)$.

Bajo este supuesto, el algoritmo define para $e$ un punto representativo $\mathbf{r}(e)$ que se calcula como el punto medio entre $\mathbf{r}(P)$ y $\mathbf{r}(Q)$, donde $P$ y $Q$ son las dos 0 -celdas que conforman la frontera de $e$. Esto asume que la geometría asociada a cada 1-celda $e$ es al menos aproximadamente la de un segmento de recta, tal como se cumple en el complejo de la figura 5.7.

Habiendo determinado $\mathbf{r}(e)$, se procede a estimar la dirección $\hat{\mathbf{v}}(e)$ del vector velocidad en este punto, usando información en $\mathcal{C}$. Un criterio simple consiste en tomar la muestra con vector posición $\mathbf{x}_{t}$ más cercano a $\mathbf{r}(e)$ y adoptar como dirección del flujo en $\mathbf{r}(e)$ al vector $\hat{\mathbf{v}}_{t}$. Alternativamente, se puede ponderar la dirección de los vectores velocidad de la nube acorde a una función peso que decrece con la distancia a $\mathbf{r}(e)$, posiblemente combinado con algún criterio para limitar el conjunto de vectores ponderados para cada $e$. Para la implementación en Python se adoptó la primera variante.

---

[comment]: # (Page break)
Habiendo determinado $\mathbf{r}(e)$ y estimado el flujo $\hat{\mathbf{v}}(e)$ en este punto, procedemos a determinar qué celda $\sigma^{*}$ de aquellas en $j(e)$ recibe el flujo a través de $e$. Para ello se adoptó el siguiente criterio:

$$
\sigma^{*}=\arg \max _{\sigma \in j(e)} \hat{\mathbf{v}}(e) \cdot \frac{\mathbf{r}(\sigma)-\mathbf{r}(e)}{\|\mathbf{r}(\sigma)-\mathbf{r}(e)\|}
$$

Dicho de otro modo, se busca maximizar el coseno del ángulo entre el vector $\hat{\mathbf{v}}(e)$ y el vector unitario que va desde $\mathbf{r}(e)$ hacia algún correspondiente $\mathbf{r}(\sigma)$. Este criterio debe producir para cada arista una única 2-celda $\sigma^{*}$. En aquellos casos en los que existan múltiples posibles candidatos para $\sigma^{*}$, en la implementación se toma solo uno de ellos, elegido arbitrariamente. El resultado es un grafo dirigido tal como se mostró en la figura 6.1, donde existe una única arista entre cada par de 2-celdas del complejo.

# 6.3. Implementación 

Un fragmento del procedimiento que construye el templex en base a una nube de puntos se muestra a continuación. El fragmento comienza aplicando el procedimiento detallado en la Sección 5 para construir una segmentación y un complejo celular en base a una nube de puntos. Luego se aplican los criterios detallados en las secciones 6.1 y 6.2 para construir el digrafo sobre el complejo. La información del digrafo queda contenida en el mapa flow_map que asigna a cada 2-celda c un conjunto de 2-celdas flow_map[c] hacia las cuales existe flujo desde c.

```
# Se construye un complejo celular en base a los datos de entrada.
cc, segmentation = OptimizingFactory(data).create_complex(kmeans_k, search_rounds, ...)
# Se obtiene la nube de puntos subyacente.
cloud = segmentation.cloud
# Se inicializa un diccionario que almacenará la información
# sobre el flujo entre caras.
flow_map = defaultdict(set)
# Se procede a evaluar todas las 1-celdas (aristas) del complejo y
# por cada una se evalúa la dirección del flujo a través de ella.
for edge in cc.edges:
    # Se obtiene la posición del punto medio de la arista,
    # asumiendo su geometría la de una recta.
    midpoint_position = edge.compute_midpoint()
    # Se estima el flujo en el punto medio
    midpoint_flow = cloud.estimate_flow_in(midpoint_position)
    # Se asegura que el flujo esté efectivamente normalizado.
    midpoint_flow = midpoint_flow / np.linalg.norm(midpoint_flow)
    # Se obtienen las caras que esta arista une.
    edge_faces = list(edge.faces)
    # Se itera a través de las caras adyacentes a la arista.
    target_face = None
    target_face_dp = -1.01
    for i, face in enumerate(edge_faces):
        # Se calcula el centroide de la cara actual.
        face_centroid = face.compute_centroid()
        # Se calcula el vector que va desde el punto medio de la
```

---

[comment]: # (Page break)
```
    # arista hacia el centroide de la cara actual.
    r = face_centroid - midpoint_position
    # Se normaliza el vector.
    r = r/np.linalg.norm(r)
    # Se calcula el producto escalar para determinar cuánto
    # del flujo es a favor de la dirección r.
    dp = np.dot(r, midpoint_flow)
    # Se mantiene registro de la cara que maximiza el flujo
    # a favor del vector r.
    if dp > target_face_dp:
        target_face, target_face_dp = face, dp
    # En el grafo de flujo se crea una arista dirigida hacia la
    # cara que maximiza la proyección del flujo sobre el vector r,
    # desde cada una de las restantes caras.
    for face in edge_faces:
        # No creamos una arista desde una cara hacia sí misma.
        if face == target_face:
            continue
        # Se mantiene registro del flujo.
        flow_map[face.order_id].add(target_face.order_id)
```

---

[comment]: # (Page break)
# Capítulo III 

## Implementación y evaluación

---

[comment]: # (Page break)
# Sección 7 

## Arquitectura de la solución

### 7.1. Introducción

Los algoritmos detallados en forma abstracta en las secciones 4,5 y 6 fueron combinados en una única herramienta de software que tiene por objetivo dar solución al problema de construir un templex en base a una nube de puntos. El código de la herramienta estará disponible tras la publicación de esta tesis en la plataforma de Git mantenida por el CIMA ${ }^{1}$. Para el desarrollo se consideraron los siguientes requisitos:

- El software está orientado a científicos e investigadores trabajando en el estudio de sistemas dinámicos no lineales. Debe ofrecer una interfaz cómoda y familiar para este segmento.
- La implementación debe ser modular tal de permitir implementar fácilmente criterios y componentes personalizados. No es necesario que el paquete sea de código cerrado.
- El software implementa el proceso de construcción de estructuras de datos sofisticadas. Debe ser posible exportar los resultados en un formato de fácil importación en entornos interactivos de desarrollo como Matlab o Mathematica.
- Es deseable poder interactuar y visualizar las estructuras de datos producidas previo a exportarlas, como sería posible en Matlab o en Mathematica. En este entorno, las estructuras de datos deben proveer interfaces de fácil utilización para la visualización de los datos.


### 7.2. Estructura de la solución

### 7.2.1. Contexto y dependencias

El lenguaje adoptado para implementar la solución es Python [28]. Este lenguaje ha ganado popularidad en el entorno científico debido a su simplicidad y alta velocidad de desarrollo. El

[^0]
[^0]:    ${ }^{1}$ https://git.cima.fcen.uba.ar/adrian.barreal/flow-templex

---

[comment]: # (Page break)
![chunk1_img-35.jpeg](chunk1_img-35.jpeg)

Figura 7.1: Diagrama de paquetes que describe el contexto en el que el código del paquete "Flow Templex" existe. Los componentes con fondo blanco son los que fueron desarrollados para este trabajo, y los componentes con fondo gris son elementos externos. Las flechas indican dependencias entre los módulos.
tipado dinámico facilita el desarrollo de algoritmos complejos, y existe una amplia variedad de paquetes orientados a la programación científica.

La estructura adoptada para el software es la de un paquete de Python de código abierto, dotado de interactividad mediante notebooks de Jupyter [29] predefinidos, que implementan distintas variantes documentadas de algoritmos que un analista podría ejecutar sobre algún conjunto de datos arbitrario que desee proveer. La comunidad científica está hoy día acostumbrada a trabajar con interfaces interactivas como la de Mathematica, y Jupyter ofrece una interfaz familiar.

Al módulo de software desarrollado se lo llamó Flow Templex y su contexto se esquematiza en la Figura 7.1. Para interactuar con el módulo, el usuario (un usuario técnico con conocimientos de programación) desarrolla sus propios notebooks de Jupyter o scripts. El paquete provee también un conjunto de notebooks demostrativos que el usuario puede usar para familiarizarse con el sistema. El paquete Flow Templex expone una interfaz programática que los programas de usuario utilizan para acceder a las funciones y a las abstracciones provistas por el módulo. El módulo a su vez utiliza paquetes estándar de Python y un conjunto de dependencias externas, las más importantes listadas a continuación:

- Numpy [30], un paquete de Python orientado al cálculo numérico. Ofrece implementaciones eficientes de vectores y matrices, así como de estructuras de mayor dimensión. También ofrece funcionalidad para operar sobre estas estructuras en forma eficiente. Las estructuras ofrecidas por Numpy son comúnmente recibidas y retornadas por otros paquetes científicos, por razones de eficiencia y compatibilidad. Por estas mismas razones es que las estructuras de Numpy fueron adoptadas como los portadores de datos fundamentales del paquete en desarrollo.
- Scikit Learn [31], un paquete que provee implementaciones de algoritmos para tareas de

---

[comment]: # (Page break)
clasificación, regresión, clustering, reducción de dimensionalidad, y otras actividades relacionadas al aprendizaje automático. Fue utilizado por su implementación de k-means.

- Networkx [32], un paquete para el trabajo con grafos. En la solución los grafos fueron implementados desde cero usando clases, pero networkx se utiliza para identificar la base mínima de ciclos al momento de construir el complejo celular. En cualquier caso, por razones de compatibilidad, los grafos específicos de la solución pueden ser convertidas a grafos de networkx para facilitar la aplicación de otros algoritmos.
- Plotly [33], un paquete para la generación de gráficos interactivos, compatible con Jupyter. Se utilizó para desarrollar las visualizaciones de las estructuras de datos desarrolladas.


# 7.2.2. Componentes internos 

El paquete Flow Templex contiene a su vez componentes que implementan la siguiente funcionalidad:

- Estructuras de datos específicas. Este componente implementa abstracciones como la nube de puntos, el complejo y el templex. Este componente construye sobre las estructuras del problema una interfaz similar o inspirada por la que networkx construye sobre el concepto de grafo. Las estructuras cuentan con correspondientes pruebas unitarias.
- Visualización. Este componente implementa herramientas compatibles con Jupyter que pueden ser usadas para visualizar las estructuras de datos específicas.
- Algoritmos utilitarios. Un conjunto de funciones utilitarias que implementan procedimientos matemáticos comunes y convenientes.
- Evaluación de resultados. Consiste en una funcionalidad que evalúa la calidad de las estructuras de datos producidas por los algoritmos, en principio en base a métricas predefinidas, pero con soporte para métricas implementadas por el usuario.
- Importación y exportación de datos. Consiste en herramientas para importar conjuntos de datos en forma apta para la utilización en la construcción de las estructuras de datos específicas. También provee herramientas para serializar y exportar las estructuras.


### 7.3. Tipos de datos e interfaz programática

En esta sección se presenta un conjunto de ejemplos que muestran cómo el paquete desarrollado simplifica la ejecución de los algoritmos detallados en las secciones anteriores. Esta simplificación se logra abstrayendo los procedimientos detrás de tipos de datos abstractos que exponen interfaces de alto nivel y encapsulan la complejidad de los procesos.

---

[comment]: # (Page break)
# 7.3.1. La nube de puntos 

El tipo de dato fundamental es la nube de puntos (tipo Cloud). El siguiente fragmento de código muestra como inicializar una nube de $n$ puntos en base a una matriz de $n \times m$ donde la $i$-ésima fila representa un vector posición asociado a la $i$-ésima muestra. El método estático from_xt asume adicionalmente que existe un ordenamiento temporal en las muestras. Bajo este supuesto se toman primeras diferencias para estimar el vector velocidad en cada punto. El tipo Cloud encapsula entonces información no solo sobre los vectores posición provistos al constructor estático, sino también los procedimientos necesarios para extrapolar y ordenar información adicional, incluyendo vectores velocidad y flujo, así como procedimientos de estimación en posiciones arbitrarias del espacio.

```
import numpy as np
from structure.point_cloud import Cloud
# Se importa un conjunto de datos.
data = np.genfromtxt('./data/rossler1.dat', delimiter=' ')
# Se inicializa un objeto tipo Cloud usando el constructor
# estático from_xt.
cloud = Cloud.from_xt(data)
# Se calcula el punto medio c entre los puntos 100 y 200
# de la secuencia.
a = cloud.position_of_point(100)
b = cloud.position_of_point(200)
c = (a+b)/2.0
# Se estima el vector flujo en el punto c.
f_c = cloud.estimate_flow_in(c)
# Se obtiene el índice del punto de la nube más cercano
# al punto c.
p_c = cloud.point_closest_to(c)
```

El tipo Cloud expone una variedad de métodos adicionales documentados en más detalle en el apéndice C.

### 7.3.2. La segmentación

Otro tipo implementado en el paquete es el tipo Segmentation, que encapsula la noción de segmentación. Este objeto internamente codifica el cómo cada punto de la nube es asignado a un segmento. Dada una nube de puntos cloud de tipo Cloud es posible construir un objeto Segmentation usando el constructor estático use_kmeans, tal como se muestra a continuación:

```
# Se segmenta la nube usando kmeans.
clustered = Segmentation.use_kmeans(cloud, k=80)
# Se obtiene el conjunto de todos los distintos colores,
# enteros que identifican segmentos, en la nube segmentada.
```

---

[comment]: # (Page break)
```
colors = clustered.colors
# Se obtienen los índices de los pseudo-centroides de cada
# segmento y se inicializa con ellos un arreglo de numpy.
pc_indices = np.array([clustered.pseudo_centroid(c) for c in colors])
# Se obtienen los vectores posición correspondientes a los
# pseudo-centroides de cada cluster.
pcs = np.array([clustered.pseudo_centroid_position(c) for c in colors])
# Se produce un gráfico 3D de la nube de puntos, asumiendo
# la nube de puntos tridimensional. Este gráfico es visible
# e interactivo en Jupyter.
clustered.draw3d(color_clusters=True)
```

Tal como lo hace el tipo Cloud, el tipo Segmentation ofrece valor agregado proveyendo funcionalidad para obtener rápida y fácilmente información sobre los distintos segmentos y sobre la nube subyacente. La clase se documenta en más detalle en el apéndice C.

# 7.3.3. El grafo de flujo 

El paquete define un tipo FlowGraph, que implementa el concepto de digrafo utilizado tanto en la construcción del complejo celular como en la construcción del simplex. El siguiente fragmento de código muestra cómo construir un digrafo dada una segmentación producida en base a una nube de puntos:

```
# Se construye la segmentación en base a una nube de puntos.
segmentation = Segmentation.use_kmeans(cloud, k=80)
# Se construye el digrafo usando el algoritmo "sequence".
fg = FlowGraph.from_segmentation(segmentation, algorithm='sequence', n=2)
```

El algoritmo sequence asume que la nube subyacente a la segmentación está ordenada cronológicamente tal de poder construir el grafo de vecindad de segmentos en base a los criterios de vecindad detallados en la Sección 5.4.

La información contenida en el digrafo depende del algoritmo de construcción usado para generarlo. El algoritmo sequence, que construye el objeto en base a una segmentación, produce un digrafo donde hay un nodo por segmento. Cada nodo tiene un identificador único igual al identificador único (el color) del correspondiente segmento.

Dado que el digrafo es un concepto relativamente genérico, el tipo FlowGraph está pensado para ser usado teniendo noción del contexto en el que fue generado. El siguiente fragmento de código muestra cómo se puede iterar por cada nodo del digrafo y calcular el centroide del correspondiente segmento:

```
# Se itera por cada uno de los nodos del grafo.
for node in fg nodes:
    # El nodo tiene un identificador único que corresponde al
    # identificador de un cluster. A continuación se obtienen
    # los vectores posición de todos los puntos en el segmento
    # correspondiente al nodo actual.
```

---

[comment]: # (Page break)
```
x = segmentation.positions_of_points_colored_as(node.id)
# Se puede obtener el centroide del segmento como se muestra
# a continuación.
x_c = segmentation.centroid(node.id)
```

El nodo del digrafo tiene también otros métodos como neighbors y neighbors_count que permiten iterar a través del correspondiente conjunto de vecinos. La clase se documenta en más detalle en el apéndice C.

# 7.3.4. El complejo celular 

El paquete provee un tipo CellComplex que implementa la noción de un 2-complejo celular. Internamente almacena un grafo donde los nodos representan 0 -celdas, y las aristas 1-celdas. El objeto almacena una construcción adicional que representa el concepto de 2-celda, a la que en este contexto se la denomina cara. El siguiente fragmento de código muestra cómo instanciar un complejo celular dado un digrafo:

```
# Se construye un complejo celular dado un digrafo que describe su estructura.
cc = CellComplex.from_flow_graph(fg)
```

El algoritmo from_flow_graph es el que se detalló en la Sección 5.5: se identifican las 2-celdas buscando una base mínima de ciclos, y los huecos se identifican buscando ciclos orientados en el digrafo.

Se define también un tipo OptimizingFactory en el mismo módulo en el que se define CellComplex. Este objeto implementa el algoritmo detallado en la Sección 5.8.

```
# Se instancia un objeto factory.
factory = OptimizingFactory(data)
# Se crea el complejo. Internamente se usa K-Means con el valor
# dado de k, y la cantidad de iteraciones de búsqueda es 10.
cc, segmentation = factory.create_complex(k=80, rounds=10)
```

La clase se documenta en más detalle en la Sección C.

### 7.3.5. El complex

El paquete implementa un tipo Ternplex que modeliza la estructura detallada en la Sección 6. El tipo Ternplex expone un constructor estático que encapsula todo el procedimiento de construcción, partiendo desde el conjunto de datos de entrada. Internamente hace uso del tipo OptimizingFactory para construir un complejo celular y sobre este generar el digrafo para producir el templex. El siguiente fragmento de código muestra cómo producir el templex dado un conjunto de datos de entrada:

```
import numpy as np
# Se carga el conjunto de datos.
data = np.genfromtxt('./data/rossler1.dat', delimiter=' ')
```

---

[comment]: # (Page break)
```
# Se construye el Templex.
    templex = Templex.from_xt(data, kmeans_k=100, search_rounds=10, quiet=False)
```

El procedimiento se muestra en el notebook de ejemplo templex.ipynb. En el apéndice C se proveen más detalles sobre el tipo Templex y su interfaz.

# 7.4. Evaluación y resultados 

### 7.4.1. Efecto del optimizador

El método create_complex del tipo OptimizingFactory recibe un parámetro booleano optativo generate_report que, de ser verdadero, le indica al algoritmo que produzca un reporte (de tipo OptimizationReport) que almacena el complejo resultante para cada ejecución, así como el correspondiente valor de aptitud.

El notebook de Jupyter optimization.ipynb muestra una variedad de pruebas que permiten observar el efecto del algoritmo de optimización. El siguiente fragmento de código es una de estas pruebas:

```
import numpy as np
from structure.cell_complex.factory import OptimizingFactory
# Se carga el conjunto de datos.
data = np.genfromtxt('./data/rossler1.dat', delimiter=' ')
cc, segmentation, report = OptimizingFactory(data).create_complex(
    k=60,
    rounds=20,
    generate_report=True,
)
scores = report.plot_scores()
```

Este fragmento construye 20 complejos (rounds=20), todos con $k=80$ y se queda con el mejor acorde a la métrica de calidad por defecto, que cuenta la cantidad de aristas necesarias para representar las líneas de juntura del complejo, y se queda con el complejo que minimiza esta cantidad. La invocación a create_complex retorna adicionalmente un reporte que guarda referencias a cada uno de los 20 complejos generados, así como al valor de aptitud correspondiente a cada uno.

La Figura 7.2 muestra un histograma donde el eje horizontal indica el negativo de la cantidad total de aristas en las líneas de juntura (la aptitud tal que a menor cantidad de aristas, mayor aptitud), y donde el eje vertical indica la cantidad de complejos producidos para cada correspondiente valor de aptitud. Se observa que, si bien la mayoría de las ejecuciones producen complejos que representan sus líneas de juntura con solo dos aristas, ocasionalmente alguna segmentación problemática produce un complejo con aptitud mucho menor. En la Figura 7.3 se muestran dos complejos: el peor complejo por un lado, y uno de los mejores complejos por otro, acorde al mencionado indicador de calidad. Usar el optimizador ayuda a mitigar el ocasional caso problemático.

---

[comment]: # (Page break)
![chunk1_img-36.jpeg](chunk1_img-36.jpeg)

Figura 7.2: Distribución de la aptitud de los complejos para una ejecución del algoritmo de optimización aleatoria con $k=80$ y 20 iteraciones. El gráfico es el que genera el método plot_scores del tipo OptimizationReport.
![chunk1_img-37.jpeg](chunk1_img-37.jpeg)

Figura 7.3: Dos de los 20 complejos producidos por una ejecución del método create_complex del tipo OptimizingFactory. La aptitud del complejo de la izquierda es -16 , siendo este el peor complejo producido por esta ejecución puntual. La aptitud del complejo de la derecha es -2, uno de los mejores complejos producidos por esta ejecución puntual. El complejo de la izquierda es evidentemente inadecuado, mientras que el complejo de la derecha es aparentemente correcto.

---

[comment]: # (Page break)
![chunk1_img-38.jpeg](chunk1_img-38.jpeg)

Figura 7.4: a) Nube de Rössler segmentada con $k=3$. b) Nube de Rössler segmentada con $k=20$. Dada la segmentación con $k=3$, el algoritmo propuesto no podrá construir un complejo homeomórfico al atractor de Rössler.

# 7.4.2. Selección del parámetro $k$ 

Bajo el supuesto de que el proceso de optimización aleatoria produce comúnmente mejores resultados que la aplicación del algoritmo una única vez, procedemos a analizar la salida producida por create_complex en función del parámetro $k$. El optimizador no trabaja determinando $k$, sino que es actualmente un parámetro que debe ser provisto por el analista. Si se provee un valor de $k$ demasiado bajo, la segmentación no será suficientemente detallada para producir un complejo adecuado. Si el valor de $k$ es demasiado alto, la distancia entre puntos sucesivos puede llegar a ser comparable al tamaño de segmentos, lo que puede a su vez producir que el algoritmo de identificación de vecinos produzca conexiones incorrectas. Idealmente, $k$ debe ser tal que los segmentos sean lo más grande posible, pero suficientemente pequeños como para que generar un vértice por segmento y conectar los vecinos produzca una representación detallada.

Consideremos por ejemplo la figura 7.4a. Esta segmentación del atractor de Rössler fue producida con $k=3$. El algoritmo detallado en la Sección 5 produce un vértice por segmento, y los conecta acorde al criterio de vecindad. Con tres vértices no es esperable producir más que un complejo homeomórfico a un triángulo, por lo que $k=3$ imposibilita producir un complejo adecuado para el atractor de Rössler.

Sin ánimo de demostración, podemos intuir en base a observación que el valor mínimo adecuado de $k$ debe ser tal de poder producir, para alguna selección inicial de centroides, una segmentación que cumpla al menos lo siguiente: en regiones visiblemente planas o similares a superficies, cualquier curva que se trace entre dos fronteras distintas del complejo debe pasar por al menos dos segmentos distintos. Consideremos por ejemplo el segmento amarillo en la Figura 7.4a. Si imaginamos al segmento como una superficie conexa, podemos trazar una curva que conecta la frontera exterior con la frontera del hueco, sin dejar el segmento. Al colapsar el segmento en un solo vértice, se pierde la dimensionalidad y la identidad de las fronteras.

Es más difícil intuir qué debe cumplir la segmentación en regiones como aquella que describe el

---

[comment]: # (Page break)
![chunk1_img-39.jpeg](chunk1_img-39.jpeg)

Figura 7.5: Complejo producido en base a la segmentación de la Figura 7.4b. El complejo no llega a capturar la línea de juntura.
![chunk1_img-40.jpeg](chunk1_img-40.jpeg)

Figura 7.6: Cantidad de aristas en líneas de juntura (eje vertical) en un complejo producido por el algoritmo de optimización aleatoria con el valor de $k$ indicado en el eje horizontal.
segmento azul de la Figura 7.4a. En este caso puntual se puede ver que, al colapsar el segmento azul a un vértice, se pierde también la reinserción de flujo que caracteriza al atractor de Rössler. Consideremos adicionalmente la segmentación de la Figura 7.4b. Esta segmentación fue producida con $k=20$. Es en apariencia más detallada que la segmentación de la figura 7.4a, pero el complejo producido (Figura 7.5) todavía no es homeomórfico al atractor de Rössler. El complejo producido en este caso es homeomórfico a una banda de Möbius, que no llega a capturar la línea de juntura.

En la figura 7.6 se muestra un gráfico de la cantidad de aristas en líneas de juntura en función de $k$, para una sucesión de ejecuciones del método create_complex de OptimizingFactory con 10 iteraciones por cada ejecución (rounds=10). El valor de $k$ aumenta en 10 en tras cada ejecución, desde $k=10$ hasta $k=200$. El experimento se encuentra también en el notebook optimization. ipynb. Observamos que la línea de juntura aparece recién con $k=30$. Una segmentación de la nube de Rössler con $k=30$ se muestra en la figura 7.7a, y el complejo producido en base a esta segmentación se muestra en la figura 7.7 b .

Para el caso puntual de este conjunto de datos, el mínimo valor de $k$ necesario para capturar consistentemente las características topológicas del atractor es aproximadamente 30 . El valor adecuado de $k$, sin embargo, depende de la nube concreta. Adicionalmente, no hay garantía de que el

---

[comment]: # (Page break)
![chunk1_img-41.jpeg](chunk1_img-41.jpeg)
(a)
![chunk1_img-42.jpeg](chunk1_img-42.jpeg)
(b)

Figura 7.7: a) Una nube de Rössler segmentada por $k$-means con $k=30$. b) Un complejo producido en base a la segmentación de la figura a.
valor de $k$ adecuado en una región sea adecuado en otra. Una mejora que podría resultar beneficioso investigar es la de aplicar una segmentación adaptativa, que utilice un valor de $k$ mayor o segmente jerárquicamente las regiones donde se detecte que esto podría ser beneficioso.

Determinar un valor de $k$ adecuado en forma sistemática para un conjunto de datos arbitrario puede ser difícil. Una estrategia para hacer frente al problema podría involucrar la determinación de un criterio adecuado de calidad que permita desarrollar un algoritmo de optimización más efectivo. Determinar este criterio de calidad y cómo evaluarlo dado una segmentación no es un problema trivial. Producir un gráfico como aquel de la figura 7.7 b puede ayudar a identificar posibles complejos candidatos.

# 7.4.3. Ejemplos de aplicación 

## El atractor de Rössler

La figura 7.8 muestra los resultados de aplicar el algoritmo a una nube de puntos generada integrando el sistema de Rössler acorde al siguiente procedimiento:

```
import numpy as np
from scipy.integrate import solve_ivp
from structure.point_cloud import Segmentation
from structure.templex import Templex
def rossler(t, r, a=0.43295, b=2.08, c=4.0):
    x, y, z = r
    dxdt = -y - z
    dydt = x + a + y
    dzdt = b + z + (x - c)
    return [dxdt, dydt, dzdt]
```

---

[comment]: # (Page break)
![chunk1_img-43.jpeg](chunk1_img-43.jpeg)

Figura 7.8: a) Nube de puntos generada muestreando una instancia del sistema de Rössler. b) Complejo celular producido en base a la segmentación de la figura 7.8a. c) Templex producido en base a la segmentación de la figura 7.8 a y el complejo $7.8 b$

```
t_span, X0 = np.linspace(0, 1000, 20000), [0.0, -2.0, 0.5]
solved = solve_ivp(rossler, [t_span[0], t_span[-1]], X0, t_eval=t_span)
data = np.vstack([solved.y[0], solved.y[1], solved.y[2]]).T
Segmentation.set_kmeans_seed(1)
templex = Templex.from_xt(data, kmeans_k=40, search_rounds=1, quiet=False)
```


# El atractor de Lorenz 

La figura 7.9 muestra los resultados de aplicar el algoritmo a una nube de puntos generada integrando el sistema de Lorenz acorde al siguiente procedimiento:

```
import numpy as np
from scipy.integrate import solve_ivp
```

---

[comment]: # (Page break)
![chunk1_img-44.jpeg](chunk1_img-44.jpeg)

Figura 7.9: a) Nube de puntos generada muestreando una instancia del sistema de Lorenz. b) Complejo celular producido en base a la segmentación de la figura 7.9a. c) Templex producido en base a la segmentación de la figura 7.9a y el complejo 7.9 b

```
from structure.point_cloud import Segmentation
from structure.templex import Templex
def lorenz(t, r, sigma=10, rho=28, beta=8/3):
    x, y, z = r
    dzdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dzdt, dydt, dzdt]
t_span, XO = np.linspace(0, 300, 20000), [0.0, -2.0, 0.5]
solved = solve_ivp(lorenz, [t_span[0], t_span[-1]], XO, t_eval=t_span).
data = np.vstack([solved.y[0], solved.y[1], solved.y[2]]).T
10 data = data[2000:]
17 Segmentation.set_kmeans_seed(1)
18
```

templex = Templex.from_xt(data, kmeans_k=70, search_rounds=1, quiet=False)

---

[comment]: # (Page break)
# 7.4.4. Dificultades adicionales con el algoritmo 

A continuación se lista una variedad de complicaciones y limitaciones identificadas durante los procesos de desarrollo y evaluación del algoritmo.

## El problema de la identificación de huecos

El algoritmo identifica huecos topológicos buscando ciclos cerrados en el grafo de flujo construido en base al algoritmo detallado en la Sección 5.6. Esto implica que aquellos huecos producidos por ramificación, como el hueco $B$ de la Figura 7.10, no serán detectados. Este es un problema fundamental con la estrategia de detección de huecos y será necesario desarrollar una estrategia alternativa si se quiere adaptar el procedimiento a casos más generales.
![chunk1_img-45.jpeg](chunk1_img-45.jpeg)

Figura 7.10: Diagrama que representa un atractor con dos huecos $A$ y $B$. En el grafo de flujo, las aristas alrededor de A conforman un ciclo cerrado y el hueco será detectado. Las aristas alrededor de B no conforman un ciclo cerrado debido a la orientación contraria de las aristas y el hueco no será detectado.

Nótese que esta no es una limitación de la teoría de grafos sino que el grafo de flujo, tal como lo genera el algoritmo detallado en esta tesis, no contiene la información necesaria para distinguir un agujero como el B en la Figura 7.10 de una 2-celda, ya que en el grafo el agujero no será necesariamente un ciclo cerrado como lo sería el agujero A. La información necesaria está (presumiblemente) en la nube de puntos y debe ser incorporada al proceso de construcción. Una posible solución puede involucrar la clasificación de los clusters generados por la segmentación distinguiendo entre aquellos que contienen agujeros y aquellos que no. Por ejemplo, bajo la hipótesis nula de que la densidad de muestras en cierto cluster es uniforme (y por lo tanto no hay hueco distinguible), ¿qué tan probable es la distribución espacial observada? Esto antes requiere contemplar el efecto que la velocidad dependiente del tiempo tiene en la densidad de muestras. Podría ser primero necesario producir una nube de puntos topológicamente equivalente a la original, que minimiza la variabilidad de la velocidad a lo largo de la trayectoria.

## El problema de la densidad no uniforme

El adecuado funcionamiento del algoritmo depende de que la nube tenga una densidad espacial de muestras aproximadamente uniforme. El problema es que en general esto depende de los parámetros del sistema. La segmentación de la Figura 7.11 fue producida integrando las ecuaciones de Rössler con parámetros $(a, b, c)=(0.43295,2.08,3.5)$. Se observa que esta nube está conformada

---

[comment]: # (Page break)
![chunk1_img-46.jpeg](chunk1_img-46.jpeg)

Figura 7.11: Nube producida integrando las ecuaciones de Rössler con $(a, b, c)=(0.43295,2.08,3.5), y$ segmentada usando $k$-means con $k=100$. Esta nube está compuesta por filamentos y dificulta la correcta ejecución del algoritmo.
aparentemente por filamentos. Por las razones detalladas en la Sección 7.4.2, el algoritmo tiene dificultades para lidiar con este tipo de nubes, produciendo resultados inevitablemente insatisfactorios.

# El problema del requisito de secuencialidad de muestras 

Existen ciertos tipos de atractores, como los atractores aleatorios, en los que no existe un ordenamiento secuencial de muestras. En estos casos no es posible aplicar los criterios de conexión de segmentos acorde al flujo. Construir complejos en este contexto más general requerirá diseñar una estrategia distinta.

### 7.4.5. Sugerencias para investigaciones posteriores

En base a la experiencia adquirida durante el diseño y la implementación del algoritmo, se proponen a continuación varias ideas que podrían llevar a procesos más robustos de construcción de complejos celulares y *templexes* en base a nubes de puntos.

## Selección de parámetros

El algoritmo propuesto requiere la definición por parte del usuario de un parámetro $k$. Si el objetivo es proveer un mecanismo robusto, general y fácil de usar, será necesario identificar un valor óptimo para este parámetro en forma automática. Implementar un procedimiento de optimización requiere definir criterios para medir la calidad de la solución. El criterio propuesto en la Sección 5.8 es uno simple, y un trabajo de investigación podría llevar a criterios más sofisticados que contemplen también la calidad de la segmentación. El criterio ideal debería permitir identificar un óptimo local cerca de aquellas soluciones que son efectivamente correctas.

---

[comment]: # (Page break)
# Generación algebraica del complejo celular 

Continuar el trabajo sobre el algoritmo que produjo el 2-complejo de la Figura 4.9 podría llevar a métodos más robustos para la construcción del complejo celular. Para generar dicho complejo se buscó producir celdas con formas acordes a las formas de los segmentos. La estrategia adoptada para lograr este objetivo consistió en identificar primero los vértices, y luego conectarlos tal de respetar siempre un conjunto de restricciones o invariantes. La implementación produjo resultados satisfactorios en regiones visiblemente planas, pero identificar, representar e implementar un conjunto adecuado de restricciones que permita lograr la correcta conexión de las celdas alrededor de las líneas de juntura es un problema que requiere más análisis. De lograr producir un algoritmo correcto, los huecos del complejo podrían ser detectados implícitamente.

## Reconstrucción del espacio topológico en base a muestras

Relacionado al problema de la selección de valores óptimos para los parámetros está el problema de encontrar una única solución dado un problema indeterminado como el de reconstruir un espacio topológico en base a muestras. Existen a priori infinitos espacios topológicos que pueden corresponder a un conjunto de muestras dado. Una forma de abordar el problema es mediante optimización, como se exploró durante esta tesis. Mediante optimización se logró una mejora notable sobre una implementación básica, como se mostró en la Sección 7.4.1. El procedimiento y la función objetivo adoptados en la implementación son primitivos, no obstante, y existe definitivamente espacio para mejora. Alternativamente, el campo de la teoría de la información ha abordado exitosamente un problema a priori similar, que es el de la reconstrucción de señales en base a muestras. Intentar aplicar técnicas de la teoría de la información a la reconstrucción del espacio topológico podría resultar productivo, aunque excede al alcance de este trabajo.

---

[comment]: # (Page break)
# Sección 8 

## Conclusiones generales

El estudio de la estructura topológica de un flujo asociado a un sistema dinámico determinístico en el espacio de estados es un tema de gran interés para la ciencia en diversos planos. Recientemente, se ha introducido un nuevo objeto matemático, denominado templex, conformado por un complejo celular y un grafo dirigido asociado a las celdas de mayor dimensión del complejo. El templex permite obtener información detallada sobre la clase topológica que condensa los mecanismos fundamentales que actúan para estructurar el flujo de un sistema dinámico. Esta tesis considera, por primera vez, el problema computacional de la implementación de la construcción algorítmica de un templex a partir de una nube de puntos.

Las metodologías previas a esta tesis consideran la construcción de complejos celulares a partir de nubes de puntos cualesquiera, sin que las coordenadas de los puntos estén necesariamente asociadas a un flujo subyacente, es decir, a la solución numérica de un sistema dinámico, o a una evolución representada en el espacio de estados a partir de series temporales. La novedad de este trabajo radica en la concepción de un algoritmo capaz de captar la organización topológica de una nube de puntos asociada a un flujo en el espacio de fases. El algoritmo desarrollado se basa en elementos de topología computacional, teoría de grafos y análisis de datos, asociando a una nube de puntos un complejo celular de celdas dotado de un digrafo que aproxima el flujo y permite extraer las propiedades topológicas del templex.

La implementación se llevó a cabo usando el lenguaje de programación Python y haciendo uso de paquetes ampliamente utilizados por la comunidad científica como numpy o scikit-learn. Se consideró que la adopción de estas tecnologías haría al producto más accesible para los potenciales usuarios interesados en esta metodología. El código fue pensado para su publicación ulterior ${ }^{1}$ con disponibilidad abierta en un formato open-source. Es importante destacar una innovación fundamental que surgió durante el desarrollo de la tesis. El algoritmo de k-means fue elegido para segmentar las nubes de puntos de modo tal que cada segmento forme una celda. La clave en este trabajo consistió en realizar la segmentación utilizando como entrada, no solamente las coordena-

[^0]
[^0]:    ${ }^{1}$ https://git.cima.fcen.uba.ar/adrian.barreal/flow-templex

---

[comment]: # (Page break)
das de los puntos de la nube, sino también las diferencias finitas en el tiempo entre puntos de las nubes, vinculadas a la existencia de un flujo subyacente. En una nube de puntos tridimensional, esto implica ejecutar k-means para seis (y no tres) coordenadas. Se trata pues de utilizar una segmentación en base a una nube de puntos $2 n$ dimensional para un espacio de fases de dimensión $n$. Hasta donde sabemos, y desde el punto de vista computacional, no hay antecedentes de este tipo de procedimiento en estudios anteriores.

Se encontró que el principal desafío en el desarrollo de un algoritmo para resolver el problema abordado en esta tesis consiste en generar variantes robustas del procedimiento, que arrojen las mismas propiedades topológicas independientemente del número de segmentos utilizados para construir el simplex y, por ende, del número y geometría de las celdas del complejo celular. Se puso especial atención en el algoritmo de construcción y pegado de los bordes de las celdas, y en la elaboración de estrategias que permitieran detectar las líneas de juntura para el flujo subyacente, denominadas joining lines. Estas líneas corresponden a las componentes de las secciones de Poincaré de un flujo y son determinantes para la correcta identificación de las propiedades topológicas de un sistema en régimen caótico. Si bien esta información puede extraerse directamente del complejo celular, es decir, sin utilizar la información provista por el digrafo, su correcta identificación es central para luego poder extraer los stripexes, es decir los caminos no redundantes que recorren el complejo celular en concordancia con el flujo.

Tal como fue implementado, el algoritmo produce resultados satisfactorios para ejemplos paradigmáticos de atractores caóticos clásicos como el de Lorenz y el de Rössler. El resultado depende, desde luego, de la elección adecuada de ciertos parámetros y de las características de bondad de la nube de puntos bajo análisis (número y distribución espacial de la colección de puntos). Si bien las estrategias adoptadas para generar los bordes de las celdas de mayor dimensión del complejo fueron diseñadas tomando como referencia nubes de puntos localmente bidimensionales en espacios de fases tridimensionales, es de esperar que el conjunto de algoritmos pueda servir de base para el desarrollo de versiones que prescindan de esta limitación. El algoritmo implementado hace uso de un procedimiento de optimización que mejora la calidad de los resultados respecto a la versión no optimizada. La optimización se funda en un criterio simple basado en invariantes topológicos que indudablemente puede ser mejorado, apelando a enfoques teóricos que quedan fuera del alcance de esta tesis. Dichos enfoques debieran apuntar a crear un conjunto de criterios adicionales para establecer el grado de confiabilidad del resultado obtenido. Por el momento, la implementación lograda requiere el criterio del especialista para determinar si los resultados producidos son satisfactorios para su ámbito de aplicación.

Esta tesis se desarrolló en el contexto de un trabajo interdisciplinario y demandó una puesta en común con aportes desde la ciencia y desde la ingeniería, cada una con una percepción distinta del problema a resolver y de los criterios de aceptación a adoptar. Desde la ingeniería informática se persiguió una solución de compromiso entre las visiones de estas disciplinas para desarrollar un marco de trabajo extensible, basado en buenas prácticas de desarrollo de software y presentado al usuario en términos de conceptos propios del dominio del problema. A su vez, el trabajo en conjunto con científicos de CIMA y CONICET bajo sus líneas habituales de investigación plantea la posibilidad de realizar publicaciones posteriores para la divulgación de los resultados aquí presentados.

---

[comment]: # (Page break)
# Bibliografía 

[1] Herbert Edelsbrunner and John L Harer. Computational topology: an introduction. American Mathematical Society, 2022.
[2] Henri Poincaré. Analysis situs. Journal de l'École Polytechnique, 1:1-121, 1895.
[3] Kolja L Kypke and William F Langford. Topological climate change. International Journal of Bifurcation and Chaos, 30(02):2030005, 2020.
[4] L.C. Kinsey. Topology of Surfaces. Lecture Notes in Computer Science. Springer-Verlag, 1993.
[5] Gisela D Charó, Christophe Letellier, and Denisse Sciamarella. Templex: A bridge between homologies and templates for chaotic attractors. Chaos: An Interdisciplinary Journal of Nonlinear Science, 32(8):083108, 2022.
[6] Robert Gilmore. Topological analysis of chaotic dynamical systems. Reviews of Modern Physics, 70(4):1455, 1998.
[7] Pablo Esteban. Científicas argentinas crean "templex": Un sistema que promete revolucionar la matemática. Página/12, Aug 2022.
[8] Denisse Sciamarella and Gabriel B. Mindlin. Topological structure of chaotic flows from human speech data. Physical Review Letters, 82(7):1450, 1999.
[9] Denisse Sciamarella and Gabriel B. Mindlin. Unveiling the topological structure of chaotic flows from data. Physical Review Letters, 64(3):036209, 2001.
[10] Gisela D Charó, Guillermo Artana, and Denisse Sciamarella. Topology of dynamical reconstructions from lagrangian data. Physica D: Nonlinear Phenomena, 405:132371, 2020.
[11] Gisela D Charó, Guillermo Artana, and Denisse Sciamarella. Topological colouring of fluid particles unravels finite-time coherent sets. Journal of Fluid Mechanics, 923:A17, 2021.
[12] Davis S. Richeson. Euler's Gem. The Polyhedron Formula and the Birth of Topology. Princeton University Press, 2008.
[13] Edward Lorenz. Deterministic Nonperiodic Flow. Journal of Atmospheric Sciences, 20(2):130148, 1963.
[14] Otto Rossler. An equation for continuous chaos. Physics Letters A, 57:397-398, 071976.

---

[comment]: # (Page break)
[15] Gisela Charó, Christophe Letellier, and Denisse Sciamarella. Templex: A bridge between homologies and templates for chaotic attractors. Chaos: An Interdisciplinary Journal of Nonlinear Science, 32:083108, 082022.
[16] Joan Birman and Robert F Williams. Knotted periodic orbits in dynamical systems i: Lorenz's equations. Topology, 22(1):47-82, 1983.
[17] Richard H Crowell and Ralph Hartzler Fox. Introduction to knot theory, volume 57. Springer Science \& Business Media, 2012.
[18] Frédéric Chazal and Bertrand Michel. An introduction to topological data analysis: fundamental and practical aspects for data scientists. Frontiers in artificial intelligence, 4:667963, 2021.
[19] Larry Wasserman. Topological data analysis. Annual Review of Statistics and Its Application, $5: 501-532,2018$.
[20] Jeff Murugan and Duncan Robertson. An introduction to topological data analysis for physicists: From lgm to frbs. arXiv preprint arXiv:1904.11044, 2019.
[21] Herbert Edelsbrunner and John Harer. Persistent homology-a survey. Contemporary mathematics, 453:257-282, 2008.
[22] Gunnar Carlsson. Topological pattern recognition for point cloud data. Acta Numerica, 23:289368, 2014.
[23] Gunnar Carlsson and Afra Zomorodian. The theory of multidimensional persistence. In Proceedings of the twenty-third annual symposium on Computational geometry, pages 184-193, 2007.
[24] Gisela D. Charó. Análisis de separación no estacionaria a partir de líneas de emisión. Facultad de Ingeniería, Universidad de Buenos Aires, 2020.
[25] Robert Ghrist. Barcodes: The persistent topology of data. BULLETIN (New Series) OF THE AMERICAN MATHEMATICAL SOCIETY, 45, 022008.
[26] Mark R Muldoon, Robert S MacKay, Jerry P Huke, and David S Broomhead. Topology from time series. Physica D: Nonlinear Phenomena, 65(1-2):1-16, 1993.
[27] Telikepalli Kavitha, Kurt Mehlhorn, Dimitrios Michail, and Katarzyna Paluch. An $\tilde{O}\left(m^{2} n\right)$ algorithm for minimum cycle basis of graphs. Algorithmica, 52:333-349, 112008.
[28] Python. https://www.python.org. Accedido: 17th Aug 2023.
[29] Jupyter. https://jupyter.org. Accedido: 17th Aug 2023.
[30] Numpy. https://numpy.org. Accedido: 17th Aug 2023.
[31] Scikit learn. https://scikit-learn.org/stable/. Accedido: 17th Aug 2023.
[32] Networkx. https://networkx.org. Accedido: 17th Aug 2023.
[33] Plotly. https://plotly.com. Accedido: 17th Aug 2023.

---

[comment]: # (Page break)
# Apéndices

---

[comment]: # (Page break)
# Apéndice A 

## Homologías

## A.1. Definiciones

Definición 7 Si $C$ es una $k$-cadena en un complejo dirigido $\mathcal{K}$, y $\partial(C)=\emptyset$, entonces $C$ es un $k$-ciclo. Al conjunto de todos los $k$-ciclos en $\mathcal{K}$ lo llamamos $Z_{k}(\mathcal{K}) \subseteq C_{k}(\mathcal{K})$.
Definición 8 Si $C$ es un $k$-ciclo en un complejo dirigido $\mathcal{K}$ tal que existe una $(k+1)$-cadena $D$ con $\partial(D)=C$, entonces $C$ es una $k$-frontera. Llamamos $B_{k}(\mathcal{K}) \subseteq C_{k}(\mathcal{K})$ al conjunto de todas las $k$-fronteras en $\mathcal{K}$.

Definición 9 Dos $k$-cadenas $C_{1}$ y $C_{2}$ en $\mathcal{K}$ son homólogas si $C_{1}-C_{2} \in B_{k}(\mathcal{K})$; es decir, si $C_{1}-C_{2}=\partial(D)$ para alguna $(k+1)$-cadena $D$. En este caso escribimos $C_{1} \sim C_{2}$.

La homología es una relación de equivalencia entre cadenas. Es reflexiva, simétrica, transitiva, y cumple que

$$
C_{1} \sim C_{2}, \text { y } C_{3} \sim C_{4} \Rightarrow C_{1}+C_{3} \sim C_{2}+C_{4}
$$

La relación de homología puede ser usada en lugar de la igualdad, dando origen a los denominados grupos de homología.
Definición 10 Sea $\mathcal{K}$ un complejo dirigido. El $k$-ésimo grupo de homología de $\mathcal{K}$ es $H_{k}(\mathcal{K})=$ $Z_{k}(\mathcal{K}) / B_{k}(\mathcal{K})$, el grupo de clases de equivalencia de elementos de $Z_{k}(\mathcal{K})$ con la relación de homología. En otras palabras, $H_{k}(\mathcal{K})$ es $Z_{k}(\mathcal{K})$ con homología utilizada en lugar de la igualdad.

## A.2. Cálculo de grupos de homología

## A.2.1. Grupos de homología del toro

En la figura A. 1 se muestra una representación simplificada de un complejo celular homeomórfico a un toro. El toro se conforma "doblando" el dibujo tal de superponer las dos instancias de la arista

---

[comment]: # (Page break)
![chunk1_img-47.jpeg](chunk1_img-47.jpeg)

Figura A.1: Un complejo celular que describe un toro.
$a$ una con otra, y las dos instancias de la arista $b$ una con otra, respetando la orientación indicada por las flechas. En este ejemplo, $\partial(\tau)=c+d-b$; por lo tanto, acorde a la definición $10, c+d \sim b$. Intuitivamente, $c+d$ y $b$ codifican información similar. Del mismo modo, $c=Q-P$ implica $Q \sim P$.

# Determinación de $H_{2}$ 

El grupo de las 2-cadenas del complejo de la figura A. 1 es generado por $\sigma$ y $\tau$. Escribimos

$$
C_{2}=\langle\sigma, \tau\rangle
$$

Una 2-cadena $C \in C_{2}$ puede ser escrita entonces como

$$
C=r_{\sigma} \sigma+r_{\tau} \tau
$$

La cadena $C$ es un ciclo si

$$
\partial(C)=r_{\sigma}(a+c+d-a-b)+r_{\tau}(c+d-b)=\emptyset
$$

Podemos re-escribir esta expresión como:

$$
\partial(C)=\left(r_{\sigma}+r_{\tau}\right) c+\left(r_{\sigma}+r_{\tau}\right) d-\left(r_{\sigma}+r_{\tau}\right) b
$$

La frontera de $C$ se anula con $r_{\sigma}=-r_{\tau}$. El conjunto de los 2-ciclos del complejo es entonces generado por $\sigma-\tau$ :

$$
Z_{2}=\langle\sigma-\tau\rangle
$$

Al no haber 3-celdas, $B_{2}=\{\emptyset\}$ y

$$
H_{2}=Z_{2} / B_{2}=Z_{2} \simeq \mathbb{Z}
$$

## Determinación de $H_{1}$

El grupo de las 1-cadenas del complejo de la figura A. 1 es generado por $a, b, c$ y $d$. Escribimos entonces

$$
C_{1}=\langle a, b, c, d\rangle
$$

---

[comment]: # (Page break)
Cualquier 1-cadena $C \in C_{1}$ puede entonces ser escrita como

$$
C=r_{a} a+r_{b} b+r_{c} c+r_{d} d
$$

Un ciclo es tal que

$$
\partial(C)=r_{a} \partial(a)+r_{b} \partial(b)+r_{c} \partial(c)+r_{d} \partial(d)=\emptyset
$$

Reemplazamos por las correspondientes fronteras de cada celda:

$$
\partial(C)=r_{a}(P-P)+r_{b}(P-P)+r_{c}(Q-P)+r_{d}(P-Q)
$$

Esto implica que $C$ es un ciclo si

$$
r_{c}(Q-P)+r_{d}(P-Q)=\left(r_{c}-r_{d}\right) Q+\left(r_{d}-r_{c}\right) P=\emptyset \Rightarrow r_{d}=r_{c}
$$

No hay condiciones impuestas sobre $r_{a}$ o $r_{b}$. Cualquier 1-ciclo del complejo es entonces un elemento del grupo

$$
Z_{1}=\langle a, b, c+d\rangle
$$

En el procedimiento de determinación de $H_{2}$ se puede observar que

$$
B_{1}=\langle c+d-b\rangle
$$

Esto implica $c+d \sim b$. Si usamos homología en lugar de igualdad, obtenemos el grupo

$$
H_{1}=Z_{1} / B_{1}=\langle a, b\rangle \simeq \mathbb{Z} \oplus \mathbb{Z}
$$

# Determinación de $H_{0}$ 

Determinamos el conjunto de las 0 -cadenas:

$$
C_{0}=\langle P, Q\rangle
$$

Por definición $\partial(P)=\partial(Q)=\emptyset$, con lo que cualquier combinación de $P$ y $Q$ es un 0 -ciclo:

$$
Z_{0}=\langle P, Q\rangle
$$

En el procedimiento de cálculo de $H_{1}$ se puede observar que

$$
B_{0}=\langle P-Q\rangle
$$

y $P \sim Q$. Con la relación de homología tenemos entonces que

$$
H_{0}=Z_{0} / B_{0}=\langle P\rangle \simeq \mathbb{Z}
$$

---

[comment]: # (Page break)
# Interpretación 

La dimensión de los grupos de homología es igual a los números de Betti $b_{k}$ del complejo. Para el toro tenemos que

- $b_{0}=\operatorname{dim}\left(H_{0}\right)=1$. Esto se interpreta como que el toro tiene una única componente conexa: es un solo cuerpo y no está compuesto por múltiples componentes aisladas.
- $b_{1}=\operatorname{dim}\left(H_{1}\right)=2$. Esto se interpreta como que el toro tiene dos ciclos planos, cada uno definido por una curva cerrada. Estos ciclos corresponden a los que se obtienen $a$ ) recorriendo el toro alrededor del eje central, y $b$ ) recorriendo el toro sobre el plano ecuatorial.
- $b_{2}=\operatorname{dim}\left(H_{2}\right)=1$. Esto implica que el toro encierra una cavidad volumétrica. Este es el volumen encerrado por la superficie de la "rosquilla."


### 1 A.2.2. Grupos de homología de la esfera

![chunk1_img-48.jpeg](chunk1_img-48.jpeg)

Figura A.2: Un complejo celular que describe una superficie esférica.
En esta sección analizamos los grupos de homología del complejo representado en la figura A.2, homeomórfico a una superficie esférica. La exposición en este caso será más breve, aunque el procedimiento de fondo es análogo al que se detalló en el cálculo de los grupos de homología del toro.

## Determinación de $H_{2}$

En este caso $\sigma$ es la única 2-celda y $\partial(\sigma)=a-a=\emptyset$. La 2-celda $\sigma$ es entonces el único generador de $Z_{2}$. Como no hay 3 -celdas tenemos que

$$
H_{2}=Z_{2}=\langle\sigma\rangle \simeq \mathbb{Z}
$$

Esto implica que el número de Betti $b_{2}=\operatorname{dim}\left(H_{2}\right)=1$; es decir, la superficie esférica encierra una cavidad volumétrica.

---

[comment]: # (Page break)
# Determinación de $H_{1}$ 

La 1-celda $a$ es el único generador de $C_{1}$. En este caso $r \partial(a)=r Q-r P$ que no se anula salvo con $r=0$. Tenemos entonces que $Z_{1}=\{\emptyset\}$. Observando la expresión $\partial(\sigma)=a-a=\emptyset$ notamos que $B_{1}=\{\emptyset\}$, por lo tanto

$$
H_{1}=Z_{1}=\{\emptyset\} \simeq 0
$$

El número de Betti $b_{1}=\operatorname{dim}\left(H_{1}\right)=0$ y la superficie esférica no tiene 1-ciclos no triviales.

## Determinación de $H_{0}$

Los generadores de $C_{0}$ son $P$ y $Q$. Cualquier combinación de $P$ y $Q$ conforma un ciclo. Notamos en $\partial(a)=Q-P$ que $Q \sim P$. El grupo $H_{0}$ es entonces

$$
H_{0}=\langle P\rangle \simeq \mathbb{Z} \quad \Rightarrow \quad b_{0}=\operatorname{dim}\left(H_{0}\right)=1
$$

## A.3. Complejos simpliciales

## A.3.1. El complejo de Čech

Definición 11 Un complejo de Čech $\check{\mathcal{C}}_{\epsilon}$, parametrizado por un número real $\epsilon$, puede ser construido en base a $\mathcal{C}$ de la siguiente forma:

1. Por cada punto $\mathbf{x} \in \mathcal{C}$, considerar una $\epsilon$-bola $\mathcal{B}_{\epsilon}(\mathbf{x})$ con centro en $\mathbf{x}$.
2. Por cada $\mathcal{S} \subseteq \mathcal{C}$, determinar $\mathcal{I}_{\mathcal{S}}=\bigcap_{\mathbf{x} \in \mathcal{S}} \mathcal{B}_{\epsilon}(\mathbf{x})$.
3. Si $\mathcal{I}_{\mathcal{S}} \neq \emptyset$, formar un $(k-1)$-símplice, con $k=|\mathcal{S}|$, con los puntos en $\mathcal{S}$.

En la Figura A. 3 se esquematiza el proceso de construcción de un complejo de Čech. Por cada punto de la nube $\{\mathrm{A}, \mathrm{B}, \mathrm{C}, \mathrm{D}, \mathrm{E}\}$ se considera una bola de radio $\epsilon$. En principio, por cada punto de la nube habrá un correspondiente 0 -símplice (un punto o vértice). La intersección entre $\mathcal{B}_{\epsilon}(\mathrm{A})$ y $\mathcal{B}_{\epsilon}(\mathrm{B})$ es no nula, por lo que se crea un 1-símplice (una arista) entre los puntos. Lo mismo ocurre con los pares $\{(\mathrm{B}, \mathrm{C}),(\mathrm{C}, \mathrm{D}),(\mathrm{D}, \mathrm{E}),(\mathrm{C}, \mathrm{E})\}$. La intersección entre $\mathcal{B}_{\epsilon}(\mathrm{C}), \mathcal{B}_{\epsilon}(\mathrm{D})$ y $\mathcal{B}_{\epsilon}(\mathrm{E})$ es también no nula, por lo que se crea un 2-símplice entre los tres vértices. El complejo resultante es el que se muestra en la Figura A.3b.

## A.3.2. El complejo de Vietoris-Rips

Definición 12 Un complejo de Vietoris-Rips $\mathcal{R}_{\epsilon}$, parametrizado por un número real $\epsilon$, puede ser construido en base a $\mathcal{C}$ de la siguiente forma:

1. Por cada $\mathcal{S} \subseteq \mathcal{C}$ determinar $d=\operatorname{máx}\left\|\mathbf{x}_{a}-\mathbf{x}_{b}\right\|$, con $\mathbf{x}_{a}, \mathbf{x}_{b} \in \mathcal{S}$.

---

[comment]: # (Page break)
![chunk1_img-49.jpeg](chunk1_img-49.jpeg)

Figura A.3: a) Una nube de cinco puntos. Sobre cada punto se dibujó una bola semi-transparente de radio $\epsilon$. b) Se construye un complejo de Čech evaluando la intersección entre las $\epsilon$-bolas.
2. Si $d \leq \epsilon$, formar un $(k-1)$-símplice con los elementos de $\mathcal{S}$, con $k=|\mathcal{S}|$.

En el complejo de Vietoris-Rips se construye un símplice por cada subconjunto de la nube que tenga diámetro máximo $\epsilon$. Entre dos puntos A y B separados por una distancia $\delta$ existirá un 1símplice si $\delta \leq \epsilon$. Entre tres puntos $\mathrm{A}, \mathrm{B}$ y C se construye un 2-símplice si no existe par de puntos separados por una distancia mayor a $\epsilon$.

---

[comment]: # (Page break)
# Apéndice B 

## Complejidad algorítmica

## B.1. Parámetros de entrada

Los siguientes son parámetros que el usuario puede especificar y son relevantes en la determinación del consumo temporal y espacial:

- $n$ : Cantidad de muestras provistas como datos de entrada. Se asume que cada muestra es un vector con una cantidad pequeña de elementos tal de poder asumir la iteración por los elementos de cada vector una operación en tiempo constante.
- $k$ : Cantidad de clusters a producir con kmeans.
- $t$ : Cantidad de iteraciones de recálculo de centroides para kmeans.
- $\tau$ : Cantidad máxima de veces a ejecutar el procedimiento de construcción del complejo celular como parte del proceso de optimización.


## B.2. Preprocesamiento de datos

La preparación de los datos tal como se implementó, incluyendo la estimación de los vectores flujo en cada punto, consiste en repetidas iteraciones lineales a través de las $n$ muestras provistas como entrada al algoritmo. El procedimiento es $\mathcal{O}(n)$ en tiempo y espacio. Los fragmentos de código más relevantes del procedimiento se listan a continuación:

```
class Cloud:
    @staticmethod
    def from_st(raw_position_sequence: np.ndarray) -> 'Cloud':
        # Take a sequence of points and extract two matrix: one that describes the spatial position
        # of the points in the cloud, and another that describes their velocity.
        position, velocity = Cloud.__extract_position_and_velocity_from_st_data(raw_position_sequence)
        # Construct a cloud from these matrices.
```

---

[comment]: # (Page break)
```
    return Cloud(position, velocity)
def __init__(self, position_data: np.ndarray, velocity_data: np.ndarray):
    self.__position_data = position_data
    self.__velocity_data = velocity_data
    self.__flow_data = Cloud.__compute_flow(velocity_data)
@staticmethod
def __extract_position_and_velocity_from_xt_data(position_data: np.ndarray):
    # We now compute velocity vector as the first-differences vector.
    velocity = np.diff(position_data, axis=0)
    # For the last position we use the last known velocity.
    velocity = np.append(velocity, [velocity[-1]], axis=0)
    # We should filter now to avoid keeping points with null velocity, since we will
    # compute flow and will be dividing by the length of the velocity.
    mask = np.linalg.norm(velocity, axis=1) > 0.0
    # Keep only the points with non-null velocity.
    return position_data[mask], velocity[mask]
@staticmethod
def __compute_flow(velocity_data: np.ndarray):
    return np.array(list(map(lambda v: v / np.linalg.norm(v), velocity_data)))
# ...
```


# B.3. Segmentación 

Segmentar requiere un procedimiento de preparación inicial que incluye la normalización de los datos y la construcción del vector de características. Este procedimiento de preparación es $\mathcal{O}(n)$ en tiempo y espacio. El elemento dominante en la etapa de segmentación es la aplicación del algoritmo kmeans. Aplicar kmeans hasta la convergencia con una cantidad variable de iteraciones es superpolinomial, pero en la implementación el procedimiento se limita a una cantidad $t$ de iteraciones tal de limitar la complejidad temporal a $\mathcal{O}(t k n)$, y la complejidad espacial a $\mathcal{O}(k+n)$. La implementación es la del paquete sklearn. El código relevante se muestra a continuación:

```
class Segmentation:
    # ...
    @staticmethod
    def use_kmeans(cloud: Cloud, k: int, max_iter=100, **kwargs):
        # Get feature matrix for segmentation. A row of the feature matrix describes one point
        # of the cloud.
        data = cloud.compute_features(**kwargs)
        # Standardize data.
        scaler = preprocessing.StandardScaler().fit(data)
        scaled = scaler.transform(data)
        # Clusterize.
        kmeans = clustering.KMeans(
            init='random',
            n_clusters=k,
            max_iter=max_iter,
            random_state=Segmentation.__kmeans_seed,
```

---

[comment]: # (Page break)
```
)
kmeans.fit(scaled)
# Get clusters.
coloring = kmeans.predict(scaled)
# Initialize a segmentation for the cloud from the data that was computed.
return Segmentation(cloud, coloring=coloring)
# ...
```


# B.4. Construcción del 1-esqueleto del complejo 

Esta etapa hacer referencia al procedimiento que construye primero un grafo de flujo entre clusters, y luego utiliza este grafo como 1-esqueleto del complejo celular. Respecto al consumo espacial, por cada color distinto se almacena un vértice en un conjunto. La cantidad máxima de vértices es igual a la cantidad de clusters, por lo que hablamos de una complejidad espacial $\mathcal{O}(k)$ para el almacenamiento de los vértices. El procedimiento incluye también la creación de un conjunto de $a$ aristas. En 2-complejos será normal observar una cantidad pequeña de vecinos por vértice, tal de poder suponer $a \propto k$. En el peor caso teórico, sin embargo, el grafo es uno completo y $a \propto k(k-1)$. Se utilizará el peor caso como cota superior aún cuando, en la práctica, un grafo completo implicaría el fracaso del algoritmo. En términos de espacio, el almacenamiento del grafo resultante es entonces $\mathcal{O}(k+a) \equiv \mathcal{O}\left(k^{2}\right)$.

Respecto al consumo temporal, el algoritmo requiere la ejecución secuencial de dos etapas diferenciadas. La primera implica un recorrido lineal visitando cada uno de los $n$ puntos de la nube. Este procedimiento es $\mathcal{O}(n)$ en tiempo. La segunda implica iterar por cada uno de los $k$ vértices del grafo que representa el 1-esqueleto del complejo, y por cada vértice visitar cada uno de sus vecinos. En el caso teórico del grafo completo, cada vértice tendrá $k-1$ vecinos y el tiempo de ejecución será $\mathcal{O}\left(k^{2}\right)$. El procedimiento conjunto es entonces $\mathcal{O}\left(n+k^{2}\right)$ en tiempo.

El fragmento de código relevante es el que se lista a continuación. El procedimiento recibe un parámetro m que puede considerarse una constante pequeña.

```
class FlowGraph!
    @staticmethod
    def from_segmentation(segmentation: Segmentation, m=2):
        # Instantiate a map that tells us, given a cluster index, to which other clusters
        # there is flow.
        flow_map = defaultdict(set)
        # Get the sequence of points.
        sequence = segmentation.cloud.position_data
        # Iterate through points in the sequence.
        for i in range(len(sequence) = m):
            # Get cluster to which this current point is assigned.
            color_i = segmentation.coloring[i]
            # Ensure that there is at least an empty set for this color.
            flow_map[color_i] = set() if color_i not in flow_map else flow_map[color_i]
            # Get clusters to which the next m points belong.
            upcoming_colors = set([segmentation.coloring[i + j] for j in range(1, m + 1)])
            # Assert that the set contains only one single color.
```

---

[comment]: # (Page break)
```
    if len(upcoming_colors) > 1:
        continue
    # Get the upcoming color from the set.
    upcoming_color_j = upcoming_colors.pop()
    # Skip if all points belong to the same cluster that the current point belongs.
    if upcoming_color_j == color_i:
        continue
    # Indicate flow from cluster i to this other cluster.
    flow_map[color_i].add(upcoming_color_j)
# Construct the flow graph from the given description.
return FlowGraph(flow_map, position=lambda c: segmentation.pseudo_centroid_position(c))
def __init__(self, flow_map: dict[int, set[int]], position: Callable):
    self.__flow_map = flow_map
    self.__position = position
    # Initialize a networkz graph from the given specification.
    self.__graph = nx.DiGraph()
    # Iterate through keys to initialize nodes.
    for node_id in flow_map.keys():
        self.__graph.add_node(node_id, position=position(node_id))
    # Iterate through key-value pairs to initialize edges.
    for node_id, neighbor_ids in flow_map.items():
        # Get position of the node with the given node_id.
        src_position = position(node_id)
        # Go through neighbors.
        for neighbor_id in neighbor_ids:
            # Get position of the neighbor node.
            dst_position = position(node_id)
            # Compute distance.
            length = np.linalg.norm(dst_position - src_position)
            # Add edge with length.
            self.__graph.add_edge(node_id, neighbor_id, length=length)
    # Initialize node wrappers. We first create a mapping from ID to node object.
    self.__nodeset = {}
    # Create the node objects without any connection.
    for node_id, node_data in self.__graph.nodes(data=True):
        # Initialize the node.
        self.__nodeset[node_id] = FlowGraphNode(nid=node_id, data=node_data)
    # Add connections.
    for node in self.nodes:
        # Get neighbors
        for neighbor_id in self.__graph.neighbors(node.id):
            node.add_neighbor(self.__nodeset[neighbor_id])
```

$\# \ldots$

# B.5. Construcción del complejo celular 

El procedimiento de construcción del complejo implica la identificación de ciclos candidatos a caras y huecos, y la distinción entre el respectivo tipo de ciclo. Este procedimiento es extenso por lo que se evitará listar código, pero consiste en la obtención de la base mínima de ciclos seguida de repetidas iteraciones a lo sumo cuadráticas sobre los objetos de dicha base, y del grafo recibido como entrada. El elemento dominante en este procedimiento es la obtención de la base mínima.

---

[comment]: # (Page break)
La implementación de networkx es $\mathcal{O}\left(a^{2} k+a k^{2} \log k\right)$ en tiempo para aristas con pesos reales, donde $a$ es la cantidad de aristas y $k$ es la cantidad de vértices [27]. Asumiendo un peor caso en el que $a \propto k(k-1)$, se tiene una complejidad temporal $\mathcal{O}\left(k^{5}\right)$. En términos de espacio, almacenar el complejo requiere almacenar una lista de listas que representará las caras del complejo. Habrá a lo sumo una cara por cada vértice del grafo del complejo tal que la cantidad de caras es $\mathcal{O}(k)$. La cantidad de vértices por cara la asumimos una constante pequeña tal que el espacio adicional requerido para almacenar las caras es también $\mathcal{O}(k)$.

# B.6. Optimización 

El procedimiento de optimización ejecutará todos los pasos anteriores una cantidad de veces $\tau$ definida por el usuario. En términos de tiempo y espacio, esto multiplica por $\tau$ el consumo temporal y espacial.

## B.7. Obtención del $\operatorname{templex}$

Para obtener el templex se construye un nuevo grafo de flujo en base al complejo y a la nube de puntos. Esto involucra una iteración por cada una de las $a$ aristas del grafo que representa al complejo. Por cada arista se visita cada uno de los $n$ puntos de la nube para encontrar el más cercano al punto medio de la arista. Esto en tiempo es $\mathcal{O}\left(k^{2} n\right)$ en el peor caso. En términos de espacio, el producto es un grafo que en consumo espacial es similar al que se utilizó para producir el complejo celular. En espacio esto es $\mathcal{O}\left(k^{2}\right)$.

El código relevante se lista a continuación:

```
class Templex:
    @staticmethod
    def from_xt(
        data: np.ndarray,
        kmeans_k=60,
        search_rounds=20,
        quiet=False):
    # Construct the cell complex from the data.
    cc, segmentation, _ = OptimizingFactory(data).create_complex(
        k=kmeans_k,
        rounds=search_rounds,
        quiet=quiet,
        generate_report=False,
    )
    # Get the underlying cloud.
    cloud = segmentation.cloud
    # To construct the templex we need a flow graph that describes the flow between
    # faces of the complex. We will first initialize an empty flow map.
    flow_map = defaultdict(set)
    # We will now evaluate all edges of the cell complex and for each one determine
    # how flow goes across it.
    for edge in cc.edges:
```

---

[comment]: # (Page break)
```
    # Get the position of the edge's midpoint.
    midpoint_position = edge.compute_midpoint()
    # Find the flow direction near the point closest to the midpoint.
    midpoint_flow = cloud.estimate_flow_in(midpoint_position)
    # Ensure that flow is actually normalized.
    midpoint_flow = midpoint_flow / np.linalg.norm(midpoint_flow)
    # Get the faces that this edge joins.
    edge_faces = list(edge.faces)
    # Go through faces attached to this edge. We will assume that only a single
    # edge can go out of each face for each edge.
    target_face = None
    target_face_dp = -1.01
    for i, face in enumerate(edge_faces):
        # Compute the vector going from the midpoint to the centroid of the face.
        face_centroid = face.compute_centroid()
        # Compute vector going from edge midpoint towards this current face.
        r = face_centroid - midpoint_position
        # Normalize.
        r = r/np.linalg.norm(r)
        # Evaluate dot product.
        dp = np.dot(r, midpoint_flow)
        # Keep track of the face towards which flow goes.
        if dp > target_face_dp:
            target_face, target_face_dp = face, dp
    # Create an edge from all other faces to the face towards flow goes.
    for face in edge_faces:
        # Strip if this is the same face.
        if face == target_face:
            continue
            # Keep track of the flow.
            flow_map[face.order_id].add(target_face.order_id)
    # New we will go through faces and keep track of the position of each one.
    face_position_data = cc.compute_face_position_data()
    # Construct now the flow graph from the given flow map.
    flow_graph = FlowGraph(flow_map, position=lambda i: face_position_data[i])
    # Return the temples object.
    return Temples(cc, flow_graph, segmentation)
def __init__(
    self,
    cell_complex: CellComplex,
    flow_graph: FlowGraph,
    segmentation: Segmentation):
    self.__cell_complex = cell_complex
    self.__flow_graph = flow_graph
    self.__segmentation = segmentation
```


# B.8. Procedimiento conjunto 

Teniendo en cuenta la secuencia conjunta de pasos, la complejidad espacial es

$$
\mathcal{O}\left(\tau\left(n+k^{2}\right)\right) \quad \text { (complejidad espacial) }
$$

---

[comment]: # (Page break)
Recolectando los elementos del consumo temporal, podemos caracterizar la complejidad temporal como

$$
\mathcal{O}\left(\tau\left(n+t k n+k^{2}+k^{5}\right)\right) \equiv \mathcal{O}\left(\tau\left(t k n+k^{5}\right)\right) \quad \text { (complejidad temporal) }
$$

---

[comment]: # (Page break)
# Apéndice C 

## Interfaz programática

En esta sección se documenta la interfaz programática de las distintas clases que hacen al paquete implementado. Todos los índices indicado en las siguientes tablas son basados en cero. Es decir, dado un índice $i$, se asume $i \geq 0$.

Cloud
El tipo Cloud representa una nube de puntos. La interfaz expone métodos para instanciar una objeto Cloud en base a datos, y para obtener luego propiedades calculadas en base a los mismos.

| Métodos estáticos |  |
| :-- | :-- |
| from_xt | Dado un arreglo de numpy de vectores posición, el método construye y <br> retorna un objeto Cloud. |
|  | @param raw_position_sequence: Arreglo de numpy que describe una se- <br> cuencia de vectores posición, orientados cronológicamente tal de poder es- <br> timar un vector velocidad tomando primeras diferencias. |
| Métodos |  |
| position_of_point | Retorna el vector posición asociado al $i$-ésimo elemento de la secuencia <br> raw_position_sequence provista como argumento al constructor estático. |
|  | @param i: Índice del vector posición a retornar, en el arreglo <br> raw_position_sequence dado a from_xt. |
| flow_at_point | Retorna el vector flujo asociado a la $i$-ésima muestra. |
|  | @param i: Índice de la muestra tal que el flujo retornado será aquel en la <br> posición position_of_point(i). |

---

[comment]: # (Page break)
| compute_features | Calcula un vector de características asociado a la $i$-ésima muestra. Este vector puede ser usado para procedimientos de segmentación, entre otros. El método puede ser redefinido por una clase heredera definida por el usuario para representar nubes de puntos en espacios de características alternativos manteniendo una interfaz homogénea, compatible con otros algoritmos del módulo. |
| :--: | :--: |
|  | @param w_pos: Flotante que define el peso asignado al vector posición. Por defecto toma el valor 1.0. |
|  | @param w_flow: : Flotante que define el peso asignado al vector posición. Por defecto toma el valor 3.0. |
| point_closest_to | Retorna el índice $i$ tal que position_of_point(i) es el más cercano posible al vector position dado como argumento. |
|  | @param position: Arreglo de numpy que especifica el vector posición del que se buscará el punto más cercano en la nube. |
| position_closest_to | Retorna position_of_point(i), donde $i$ es valor retornado por point_closest_to con argumento igual a position. |
|  | @param position: Arreglo de numpy que especifica el vector posición del que se buscará el punto más cercano en la nube. |
| estimate_flow_in | Estima el vector flujo en la posición position dada como argumento. El algoritmo por defecto retorna el vector flujo correspondiente a la $i$-ésima muestra, donde $i$ es el valor retornado por point_closest_to(position). El método puede ser redefinido por una clase heredera para implementar formas alternativas de estimar el flujo en una posición dada, en base a información contenida en la nube. |
|  | @param position: Arreglo de numpy que especifica un vector posición donde estimar el flujo. |
| Propiedades |  |
| position_data | Retorna el arreglo de numpy provisto como argumento a from_xt, tal que la $i$-ésima fila es el vector posición de la $i$-ésima muestra, igual a position_of_point(i). |
| velocity_data | Retorna el arreglo de numpy donde la $i$-ésima fila es el vector velocidad asociado a la $i$-ésima muestra. |
| flow_data | Retorna el arreglo de numpy donde la $i$-ésima fila es el vector flujo asociado a la $i$-ésima muestra. |

---

[comment]: # (Page break)
| Segmentation |  |
| :--: | :--: |
| Encapsula el procedimiento de segmentación de una nube de tipo Cloud, y los medios para interactuar con la información asociada. |  |
| Métodos estáticos |  |
| set_kmeans_seed | Define la semilla a utilizar para la segmentación con k-means. Esto hace a la segmentación determinística. |
|  | @param seed: Entero. La semilla a establecer. |
| use_kmeans | Método estático que produce una segmentación en k clusters usando k-means, dada una nube de puntos Cloud, con un máximo de iter iteraciones. |
|  | @param c: Objeto tipo Cloud. La nube de puntos. |
|  | @param k: Entero. El valor de $k$ a usar con k-means. |
|  | @param iter: Entero. Cantidad máxima de iteraciones. Es un parámetro de k-means. Por defecto toma el valor 100. |
|  | @param *kwargs: Argumentos adicionales. Serán provistos a compute_features para obtener el vector de características de la nube. |
| from_centroids | Método estático que inicializa una segmentación usando k-means, con un conjunto de centroides predefinido provisto como una colección de índices que hacen referencia a muestras en la nube provista como primer argumento. |
|  | @param c: La nube de puntos a segmentar. |
|  | @param centroids: Un vector de enteros donde cada entero es el índice de una muestra en la nube. La posición de la correspondiente muestra estará en el conjunto inicial de centroides. |
| Métodos |  |
| points_colored_as | Retorna el conjunto de los índices de puntos en la nube que han sido ubicados en el cluster cuyo identificador único es color. |
|  | @param color: Entero que identifica al cluster cuyos índices de muestras serán retornados. |

---

[comment]: # (Page break)
| positions_of_points_ <br> colored_as | Análogo al método points_colored_as, excepto que esta variante retorna <br> un arreglo de numpy con los vectores posición de cada muestra en el cluster <br> cuyo identificador único es color. |
| :-- | :-- |
|  | @param color: Entero que identifica al cluster cuyos vectores posición <br> serán retornados. |
| pseudo_centroid | Dado un conjunto de vectores posición correspondientes a muestras en un <br> cierto cluster, el centroide de los mismos es el primer momento del conjunto. <br> El centroide no es necesariamente una muestra de la nube, sin embargo, y <br> puede no ser útil en casos en los que necesitamos un vector representativo <br> que sí corresponda a una muestra. Llamamos entonces pseudo-centroide al <br> punto de la nube más cercano al centroide geométrico. Este método retorna <br> el índice en la nube del pseudo-centroide del cluster cuyo identificador único <br> es el entero color. |
| pseudo_centroid_position | @param color: Entero que identifica al cluster del que se retornará el <br> índice en la nube de su pseudo-centroide. |
|  | Análogo al método pseudo_centroid, con la variante de que este método <br> retorna el vector posición de la muestra y no el índice en la nube. |
|  | @param color: Entero que identifica al cluster del que se retornará el <br> vector posición del pseudo-centroide. |
|  | Retorna el vector posición del centroide geométrico correspondiente al clus- <br> ter con identificador único color. |
| draw3d | @param color: Entero que identifica al cluster del que se retornará el <br> vector posición correspondiente al centroide. |
|  | Método utilitario que delega al módulo de visualización para producir un <br> gráfico visible en notebooks de Jupyter. |
|  | @param show: Conjunto de enteros. El conjunto de clusters que serán <br> representados en el gráfico. Si el argumento es nulo se mostrarán todos los <br> clusters. |
|  | @param color_clusters: Es un booleano tal que, si el valor es falso, todos <br> los clusters serán representados con el mismo color; de ser verdadero, a cada <br> cluster se le asignará un color de entre un conjunto de colores por defecto. |

---

[comment]: # (Page break)
| cloud | Retorna la nube de puntos en la cuál se basó la segmentación, la misma <br> que debe haber sido provista como primer argumento a use_kmeans <br> o from_centroids. |
| :-- | :-- |
| coloring | Retorna una lista donde el $i$-ésimo elemento es un identificador único <br> del cluster al que pertenece la $i$-ésima muestra de la nube (aquella <br> cuya posición es cloud.position_data(i)). |
| colors | Retorna el conjunto de todos los identificadores únicos de clusters <br> producidos por la segmentación. Es decir, es el conjunto de todos los <br> elementos distintos en la lista retornada por el método coloring. |
| color_count | Retorna la cantidad de colores distintos. Es igual a la cantidad de <br> elementos en el conjunto retornado por el método colors. |


| FlowGraph |  |
| :-- | :-- |
| Implementa los procedimientos necesarios para construir y trabajar con grafos de flujo. Se trata <br> sencillamente de un grafo dirigido, con métodos estáticos utilitarios que permiten construirlo en <br> base a otros objetos de dominio, incluyendo objetos tipo Segmentation. |  |
| Métodos estáticos |  |
| from_segmentation | Construye y retorna un objeto FlowGraph dada una segmentación s. |
|  | @param s: Objeto tipo Segmentation. La segmentación en base a la cual <br> construir el grafo de flujo. El algoritmo asume que las muestras de la nube <br> s. cloud están ordenadas cronológicamente tal que es posible estimar la <br> dirección del flujo observando muestras sucesivas. |
|  | @param n: Entero tal que dos segmentos $\mathcal{A}$ y $\mathcal{B}$ serán conectados por una <br> arista dirigida del primero hacia el segundo si existe alguna muestra $s$ en <br> $\mathcal{A}$ tal que las siguientes $n$ muestras a $s$ estén todas en $\mathcal{B}$. |
| Métodos | Convierte al objeto FlowGraph en un grafo dirigido de networkx de <br> tipo nx.DiGraph. |
| to_networkx | Retorna un vector posición asociado al $i$-ésimo nodo, acorde a un criterio <br> dependiente de cómo se inicializó el grafo. |
| get_position |  |

---

[comment]: # (Page break)
|  | @param i: Índice del nodo en el grafo de flujo. Si se usó el método estático from_segmentation(s), la posición del nodo i corresponde al pseudocentroide del i-ésimo segmento, el que se obtiene invocando al método s.pseudo_centroid_position(i). |
| :--: | :--: |
| draw3d | Un método proxy que simplifica el acceso al módulo de visualización. Los parámetros adicionales kw serán provistos al método draw3d del módulo de visualización. <br> @param position: Función que recibe un único parámetro entero. Es un criterio que, dado un índice i de un nodo del grafo, retorna un vector posición. Si no se provee ninguno, el valor por defecto es análogo a invocar get_position(i). |
| Propiedades |  |
| flow_map | Retorna un diccionario d tal que, dado un identificador único c de un cluster $\mathcal{A}$ (es decir, el "color" del cluster) referencia un conjunto d[c] que es el conjunto de todos los segmentos hacia los cuales existe flujo desde $\mathcal{A}$. |
| node_count | Retorna la cantidad de nodos en el grafo de flujo. |
| nodes | Retorna un iterador a través de los nodos del grafo de flujo. |


| CellComplex |  |
| :-- | :-- |
| Implementa los procedimientos necesarios para construir y trabajar con complejos celulares orien- <br> tados. |  |
| Métodos estáticos | Método estático que construye un complejo celular en base a un grafo <br> de flujo. El algoritmo calcula la base mínima de ciclos en términos de la <br> longitud total del ciclo, y aplica los criterios detallados en la Sección 5 para <br> distinguir caras de huecos y luego orientarlas. Si los nodos del grafo de flujo <br> tienen información geométrica asociada, las celdas del complejo también <br> la tendrán. |
|  | @param g: Objeto tipo FlowGraph. El grafo de flujo en base al cual cons- <br> truir el complejo. |

---

[comment]: # (Page break)
| to_networkx | Retorna un par de elementos: una representación del complejo en for- <br> ma de grafo de networkx, y un objeto Faces que no es más que un alias <br> de un tipo list[1ist[int]], y representa una lista de las "caras" o <br> 2-celdas del grafo identificadas por el algoritmo de construcción. Cada <br> 2-celda se representa como una lista de identificadores de nodos en el <br> grafo, y recorrer la lista correspondiente a una 2-celda es análogo a <br> un recorrido sobre la frontera de la 2-celda. |
| :-- | :-- |
| compute_face_position_data | Calcula y retorna una matriz de numpy donde la $i$-ésima fila es un <br> vector posición representativo de la $i$-ésima cara, la que se obtiene <br> iterando hasta el índice $i$ sobre el iterador resultante de invocar al <br> método faces. El vector posición representativo de una cara es sim- <br> plemente el centroide de los vectores posición representativos de los <br> nodos que conforman su frontera. |
| draw3d | Método que dibuja el complejo celular delegando al módulo de visu- <br> lización. |
| Propiedades |  |
| joining_line_edge_count | Retorna la cantidad de aristas usadas para representar líneas de jun- <br> tura. En otras palabras, es la cantidad total de aristas que conforman <br> la frontera de tres o más 2-celdas. |
| nodes | Retorna un iterador sobre los nodos del grafo que define el 0-esqueleto <br> del complejo. Los nodos son de tipo Node y su interfaz se describirá <br> luego. |
| edges | Retorna un iterador sobre las aristas del grafo que define el 1-esqueleto <br> del complejo. Las aristas son de tipo Edge y su interfaz se describirá <br> luego. |
| faces | Retorna un iterador sobre las "caras" del grafo que define el 2- <br> esqueleto del complejo. Las caras son de tipo Face y su interfaz se <br> describirá luego. |


| Node |  |
| :-- | :-- |
| Representa un nodo de un grafo. Es inicializado por el algoritmo de construcción del complejo <br> celular. |  |
| Propiedades |  |
| id | Entero que identifica unívocamente al nodo en el grafo. |

---

[comment]: # (Page break)
| position | Arreglo de numpy. Vector posición asociado al nodo. |
| :-- | :-- |
| neighbors | Iterador sobre todos los nodos vecinos. |


| Edge |  |
| :-- | :-- |
| Representa una arista en un grafo. Es inicializada por el procedimiento de construcción del complejo <br> celular. |  |
| Métodos |  |
| compute_midpoint | Retorna un arreglo de numpy que representa el vector posición del pun- <br> to medio entre las posiciones de los correspondientes nodos que la arista <br> conecta. |
| Propiedades |  |
| id | Cadena de caracteres. Es un identificador único generado en base a <br> los identificadores únicos de los correspondientes nodos que esta arista <br> conecta. |
| nodes | Retorna una tupla de dos objetos de tipo Node, los nodos conectados <br> por esta arista. |
| faces | Retorna un iterador sobre una colección de objetos tipo Face. Esta es <br> la colección de todas las caras o 2-celdas del complejo que en su fron- <br> tera tienen a la arista en cuestión. Para aristas en líneas de juntura, <br> la cantidad de resultados será mayor a 2. |
| face_count | Entero. La cantidad de cadas o 2-celdas del complejo que en su fron- <br> tera tienen a la arista en cuestión. Para aristas en líneas de juntura, <br> el valor de esta propiedad es mayor a 2. |


| Face |  |
| :-- | :-- |
| Representa una cara o 2-celda de un 2-complejo. Es inicializada por el procedimiento de construcción <br> del complejo celular. Internamente se representa como una lista ordenada de nodos (Node). El orden <br> de la lista define un recorrido que le da orientación a la cara. |  |
| Métodos |  |
| to_list | Retorna una lista de enteros, donde el $i$-ésimo entero es el identificador <br> n.id del $i$-ésimo nodo n. |

---

[comment]: # (Page break)
| compute_centroid | Arreglo de numpy. Es el vector posición calculado como el centroide de todos los nodos que hacen a la frontera de la cara. |
| :--: | :--: |
| has_edge | Booleano. Verdadero si en el camino que define la frontera de la cara existe una arista entre los nodos a y b. |
|  | @param a: Objeto tipo Node. |
|  | @param b: Objeto tipo Node. |
| Propiedades |  |
| id | Identificador único de la cara, generado en base a los identificadores únicos de los nodos que componen la frontera de la cara. Es independiente del orden de los nodos. |
| order_id | Un identificador único de la cara provisto por el algoritmo de construcción del complejo al momento de instanciarla. |
| nodes | Retorna un iterador sobre una colección de objetos Node, los nodos que conforman la frontera de la 2-celda. |
| node_ids | Similar a nodes, salvo que el iterador es en este caso sobre la secuencia de identificadores de nodos (node.id) en vez de tratarse de los objetos Node mismos. |
| node_count | Entero. Retorna la cantidad de nodos que hacen a la frontera de la cara. |


| Templex |  |
| :--: | :--: |
| Métodos estáticos |  |
| from_xt | Construye un Templex dado un conjunto de datos data. Internamente utiliza el tipo OptimizingFactory para producir primero un complejo mediante optimización aleatoria (Sección 5.8) y para luego construir el Templex acorde al procedimiento detallado en la Sección 6. |
|  | @param data: El conjunto de datos, como arreglo de numpy. Debe ser apto para ser provisto al método from_xt del tipo Cloud. |
|  | @param kmeans_k: Entero. El valor de $k$ a proveer para k-means. |

---

[comment]: # (Page break)
|  | @param search_rounds: Entero. La cantidad de rondas para el algoritmo <br> de optimización aleatoria. |
| :-- | :-- |
|  | @param quiet: Booleano. De ser falso, se emitirán mensajes al inicio de <br> cada ronda. |
| Métodos |  |
| to_json | Serializa el Templex a un archivo JSON. |
| draw3d | Un método proxy que simplifica el acceso al módulo de visualización. Pro- <br> duce una representación gráfica del templex, visible en la interfaz de Jupy- <br> ter. |
| Propiedades |  |
| id | Identificador único de la cara, generado en base a los identificadores <br> únicos de los nodos que componen la frontera de la cara. Es indepen- <br> diente del orden de los nodos. |
| cloud | Objeto de tipo Cloud. La nube de puntos producida. |
| segmentation | Objeto de tipo Segmentation. La segmentación producida y utilizada <br> durante el proceso de construcción del complejo celular. |
| cell_complex | Objeto de tipo CellComplex. El complejo celular. |
| flow_graph | Objeto de tipo FlowGraph. El grafo que describe el flujo entre celdas <br> del complejo celular. |

---

[comment]: # (Page break)
