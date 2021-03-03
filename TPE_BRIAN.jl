### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 09062294-3e5f-11eb-176f-dfcbf841f111
begin
	import Pkg
	Pkg.add.(["DSP", "FFTW", "FileIO", "WAV", "PlutoUI", "Plotly"])
	Pkg.add.(["LibSndFile", "SampledSignals", "Statistics"])
	Pkg.add(["ImageMagick", "Images", "ImageIO", "Plots", "IterTools"])
	#Pkg.add("Calculus")
	using Plots
		plotly()
	using DSP, FFTW
	using DSP.Windows: hamming, rect, hanning
	using PlutoUI
	using WAV
	using FileIO
	using Images
	using LibSndFile
	using Statistics
	using SampledSignals
	
	
	"Pad vector with zeros on the right until its length is `n`"
	padright(x, n) = copyto!(zeros(eltype(x), n), x)
	
	"Plot discrete functions"
	stem(args...; kwargs...) = sticks(args...;
									marker=:circle,
									leg=false,
									kwargs...)

	stem!(args...; kwargs...) = sticks!(args...;
									marker=:circle,
									leg=false,
									kwargs...)

	"""
	Función módulo pero con offset (opcional)
	Manda a `t` al intervalo [from, from+length)
	sumándole o restándole algún múltiplo de `len`
	"""
	cshift(t, len, from=0) = mod(t - from, len) + from
	

	# Espectrograma
	using IterTools
	function stft(x; overlap, window, nfft, rest...)
	  nwin = length(window)
	  @assert overlap < nwin

	  res = [ fft(padright(xseg .* window, nfft))
		for xseg in partition(x, nwin, nwin - overlap)]

	  return [ res[i][j] for j in 1:nfft, i in eachindex(res)]
	end


	specplot(x::AbstractMatrix; kwargs...) = 
		@error "You are entering a Matrix (2D Array). I need a Vector (1D Array)."
	function specplot(x::AbstractVector;
		  fs=1,
		  onesided=false,
		  xaxis="Tiempo (s)",
		  yaxis="Frecuencia (Hz)",
		  window=hamming(div(length(x), 16)),
		  overlap=0.5,
		  nfft=length(window),
		  kws...)

		window isa Integer && (window = rect(window))
		overlap isa AbstractFloat && (overlap = round(Int, length(window) * overlap))

		mat = stft(x; overlap=overlap, window=window, nfft=nfft)

		fmax = fs
		if onesided
		  mat = mat[1:div(size(mat, 1) + 2, 2), :]
		  fmax = fs/2
		end

	  toffset = length(window) / 2fs
	  times = range(toffset; length=size(mat, 2), stop=length(x)/fs - toffset)
	  freqs = range(0; length=size(mat, 1), stop=fmax)

		# Reubico las frecuencias negativas arriba de todo
	  if !onesided
		freqs = cshift.(freqs, fs, -fs/2)
		ord   = sortperm(freqs)
		mat   = mat[ord, :]
		freqs = freqs[ord]
	  end

		return heatmap(times, freqs, log.(abs.(mat) .+ eps());
			  xaxis=xaxis, yaxis=yaxis,
			  seriescolor=:bluesreds, legend=true, kws...)
	 return times, freqs, mat
	end

	function specplot(x :: AbstractVector{<:AbstractFloat}; kws...)
		return specplot(convert.(Complex, x); onesided=true, kws...)
	end
	
	# Polos y ceros
	
	zeropolegain(pr) = DSP.ZeroPoleGain(pr)
	zeropolegain(z, p, g) = DSP.ZeroPoleGain(z, p, g)
	polynomialratio(zpg) = DSP.PolynomialRatio(zpg)
	function polynomialratio(b, a)
	  n = max(length(a), length(b))
	  return DSP.PolynomialRatio(padright(b, n), padright(a, n))
	end
	getpoles(zpg) = DSP.ZeroPoleGain(zpg).p
	getzeros(zpg) = DSP.ZeroPoleGain(zpg).z
	getgain(zpg) = DSP.ZeroPoleGain(zpg).k
	getnumcoefs(pr) = trimlastzeros!(reverse(DSP.PolynomialRatio(pr).b.coeffs))
	getdencoefs(pr) = trimlastzeros!(reverse(DSP.PolynomialRatio(pr).a.coeffs))
	function trimlastzeros!(a)
	  !iszero(a[end]) && return a
	  pop!(a)
	  return trimlastzeros!(a)
	end



	function zplane(zs, ps; kwargs...)
		scatter(real.(zs), imag.(zs);
			  marker = (:black, :circle), label="Cero", kwargs...)
		scatter!( real.(ps), imag.(ps);
			marker = (:red, :xcross), label="Polo", kwargs...)
	  ts = range(0,stop=2pi;length=100)
	  plot!(cos.(ts), sin.(ts); aspect_ratio = 1, legend=false, kwargs...)
	end

	zplane(pr::DSP.PolynomialRatio; kwargs...) = 
		zplane(DSP.ZeroPoleGain(pr); kwargs...)
	
	DSP.filt(zpg::DSP.ZeroPoleGain, r...; kwargs...) = 
		filt(polynomialratio(zpg), r...; kwargs...)



	# Delta
	d(n) = n == 0 ? 1. : 0.

	# Escalón
	u(n) = n >= 0 ? 1. : 0.


end;

# ╔═╡ 8f1394ee-3e63-11eb-0093-e75468460dc5
md"# Trabajo Práctico Especial de Señales y Sistemas"

# ╔═╡ 7c04611e-3e61-11eb-0aa5-eb97132ace53
html"""
<h2 class="pm-node nj-subtitle">Identificación de canciones mediante huellas digitales acústicas</h2>
"""

# ╔═╡ 82bdb9f0-7a26-11eb-0df3-0553dd890af1
md"""
- **Alumno: Brian Alex Fuentes Acuña.**

- **Padron: 101785.**

- **Curso: Rui Rojo.**
"""

# ╔═╡ adc46380-3e63-11eb-2422-5bfe1b5052ba
md"""
# Introducción

Imaginen la siguiente situación: se encuentran en un bar, tomando su trago favorito, inmersos en esa interesante conversación cuando de repente . . .  “pará, pará, ¿qué es eso que suena?”. Entre el ruido general alcanzan a distinguir esa canción que no escuchaban en tanto tiempo, que habían obsesivamente intentado encontrar pero que, a pesar de ello, nunca llegaron a dar siquiera con el nombre del intérprete. Su corazón se estremece . . .  Ahora que la tienen ahí, sonando nuevamente, la situación es bien diferente a la de hace algunos años: toman su teléfono celular, abren alguna de las aplicaciones de reconocimiento de audio que tienen y, en cuestión de segundos, aparece la información completa del tema, la banda, ¡e incluso se atreve a sugerir artistas similares!

Si el protagonista de esta ficción fuera estudiante de Señales y Sistemas, no debería extrañarnos que surjan en él varios interrogantes: ¿cómo hizo el programa para reconocer la canción escuchando sólo 10 segundos?, ¿cómo puede funcionar con el ruido de un bar de fondo? y ¿cómo pudee ser que además identifique tan rápido a seta banda que nadie conoce?

Resulta que todo este proceso se basa en reconocer *huellas digitales acústicas*, una técnica robusta y eficiente que puede entenderse en términos de Señales y Sistemas.

## Reconocimiento mediante huellas digitales acústicas

El reconocimiento de audio mediante huellas digitales acústicas es una técnica que busca identificar una pieza de audio, contrastando la información contenida en dicha pieza contra la almacenada en una base de datos de pistas conocidas. Esta técnica comienza a desarrollarse desde el año 2000, pero es durante la última década que tomé mayor impulso, siendo posible hoy en día encontrar una variedad de aplicaciones para smartphones capaces de identificar casi cualquier pieza de audio en segundos \[4\].

Existen varias formas dea bordar el problema del reconocimiento, pero todas persiguen la misma finalidad:

* Simplicidad computacional: el reconocimiento debe realizarse en forma rápida.
* Eficiencia en el uso de memoria y buena capacidad de discriminación: existen aproximadamente 100 millones de canciones grabadas \[1\].
* Robustez ante degradaciones: ruido de fondo, degradación por compresión digital del audio, ecualización por respuesta en frecuencia no plana del lugar y/o parlantes, etc.
* Granularidad: capacidad de reconocimiento utilizando sólo un segmento del audio.
* Alta tasa de aciertos, baja tasa de falsos positivos y de rechazos.

El procedimiento básico consiste en analizar el segmento de audio buscando extraer características particulares en el esapcio tiempo-frecuencia (es decir, en su espectrograma) que sirvan para identificarlo luego. Una característica podría ser, por ejemplo, la potencia que posean las diferentes bandas de frecuencias. Al conjunto de estas características se las denomina *huella digital acústica*.

![sys-tpe1.png](https://i.imgur.com/FoKz8Tw.png)

Este procedimiento de extracción de huellas se utiliza tanto para confeccionar la base de datos de canciones conocidas como para posteriormente reconocer segmentos de audio que se desean identificar consultando la base de datos. 

La elección de las características es la base del funcionamiento del método. Deben ser lo suficientemente específicas como para identificar a cada canción unívocamente pero a la vez ocupar poca memoria para permitir realizar la búsqueda en forma rápida y eficiente. También deberán ser relativamente inmunes ante ruidos y distorsiones del audio de manera de lograr un sistema de reconocimiento robusto. 

Diferentes características han sido propuestas con el fin de lograr estas especificaciones. Algunos ejemplos que se extraen del dominio tiempo-frecuencia son los MFCC \[4\] (Mel-Frequency Cepstrum Coefficients, comúnmente utilizados para la representación del habla), SFM \[1\] (Spectral Flatness Measure, una estimación de cuánto se aproxima una banda espectral a ser un tono o ruido), la ubicación de los picos de potencia del espectrograma (la base del algoritmo de Shazam \[9, 8\]), la energía de las bandas \[6\], etc. Además del espectrograma, otras técnicas se han utilizados como *wavelets* y algoritmos de visión \[2\], algoritmos de machine learning \[3\], y meta-características \[7\], útiles en grabaciones muy distorsionadas. 

En este trabajo desarrollaremos una implementación simple que utiliza como características el signo de la diferencia de energía entre bandas \[5\]. 

Cabe destacar que esta característica sólo sirve para reconocer las mismas versiones de las canciones que se almacenan en la base de datos: es decir, no está diseñada para reconocer versiones en vivo o interpretaciones diferentes a la original.

# Descripción del sistema de reconocimiento

El sistema de reconocimiento consta de dos bloques principales:

1. El algoritmo para extraer las huellas digitales acústicas.
2. La base de datos que almacena las huellas acústicas de las canciones conocidas y permite realizar búsquedas a partir de la huella acústica de una canción desconocida.

El algoritmo para extraer las huellas acústicas toma como entrada el archivo con el audio, lo pre-procesa, extrae las características, y entrega como resultado la huella acústica.


![sys-tpe2.png](https://cdn.nextjournal.com/data/QmRYiHJKoMuXyNybtpYDvExYoxBiwoA4dHQei9Quar5Sit?filename=sys-tpe2.png&content-type=image/png)

La base de datos guardará las huellas acústicas de las canciones conocidas. Tiene además un sistema de búsqueda (en inglés query) tal que al realizar un query con una huella acústica dada nos devuelve la canción – más probable – a la cual correseponde.

El esquema del sistema completo se presenta a continuación:

![sys-tpe3.png](https://cdn.nextjournal.com/data/QmYQPVLaSUHo9XBXL2m3HqMnbA2QAY9jWg7VqoFeby8N77?filename=sys-tpe3.png&content-type=image/png)

Observe que el algoritmo de extracción de huellas se utiliza tanto en la creación de la base de datos de canciones conocidas como para el reconocimiento de audios desconocidos.

## Algoritmo de extracción de huellas digitales acústicas

El algoritmo de extracción de huellas acústicas tiene como finalidad extraer características a partir del espectrograma del audio. Primeramente, se acondiciona el audio con un pre-procesamiento. Luego se realiza un espectrograma del audio pre-procesado y finalmente se extraen las características para cada intervalo de tiempo.

### 1. Pre-procesamiento del audio

Se comienza convirtiendo el audio, que suele ser estar en estéreo, a un audio monocanal, promediando los canales. 

Luego, se reduce su frecuencia de muestreo, aprovechando que son las frecuencias bajas las que contienen la información más útil para la extracción de características. En general, el audio está muestreado a 44100 Hz (calidad de CD), pero para este algoritmo de reconocimiento basta con tenerlo muestreado a 1/8 de su frecuencia original, 5512.5 Hz. Esto permite trabajar con menos muestras, aliviando la carga computacional. Para realizar los siguientes ejercicios, utilice el archivo `Pink.ogg`.

#### Ejercicio 1)

**Cargue la pista de audio. Verifique que la variable cargada es una matriz con 2 columnas, cada una correspondiendo a un canal de audio. Promedie los canales utilizando al función mean para tener 1 solo canal y grafique una porción de la señal resultante. Escuche el audio resultante para verificar el resultado.**
"""

# ╔═╡ a3bf22c4-3ea3-11eb-3d3d-adfdfc171c33
# La frecuencia de muestreo está fija a 44100 Hz
const sr = 44100

# ╔═╡ 54e21790-7a14-11eb-14f3-9bf2cbf17fe6
md""" Se carga el audio a procesar y se obtienen sus muestras en 2 columnas, canales estéreo."""

# ╔═╡ d132a762-3ea3-11eb-3494-692576a31f34
function load_audio(fn)
	
    xraw = load(fn)
	
		# Chequeo que la frecuencia de muestreo sea `sr`
    @assert samplerate(xraw) == sr  samplerate(xraw) 

    return collect(xraw)
end

# ╔═╡ aeea5400-794c-11eb-35b4-b763196181af
f, col_Pink= size(load_audio("Pink.ogg"))

# ╔═╡ 4612d502-7961-11eb-0a80-11d812178ad8
md""" Se observa que cargando el audio se obtienen 441000 muestras para dos canales, sus columnas.

Para obtener un solo canal se crea la función **to_mono()** quien recibe el audio con sus dos canales retornando un vector con las muestras promediadas
."""

# ╔═╡ 28c5ed26-3e6b-11eb-1d44-01e209b92f00
#Recibe las muestras del vector y las convierte en un audio para asi poder oirlo.
sound(x) = SampleBuf(x, sr);

# ╔═╡ 92649f90-73cd-11eb-0df8-4958d753607d
#toma un audio (.ogg) de dos canales, los promedia a un solo canal.
function to_mono(fn)
	stereo = load_audio(fn)
	mono = vec( mean(stereo, dims=2) )
	return mono
end


# ╔═╡ a2fa88b0-73cd-11eb-1336-9fbf72b0ddd8
tst = to_mono("Pink.ogg");


# ╔═╡ b65f8000-7962-11eb-1ba6-213e0aac1846
md""" Se puede escuchar el audio mono canal obtenido haciendo uso de la función **sound()** defina previamente."""

# ╔═╡ a7727c3e-73cd-11eb-3bdc-8dd64b6b43ad
sound(tst)


# ╔═╡ a83a1200-73cd-11eb-353c-751703d316cd
begin
  to_samp = 1*sr;     #quiero graficar el prmer segundo.
  n=1:1:to_samp;
  tseg = n./sr
end;

# ╔═╡ a56414d0-7621-11eb-0140-e57ff8c679ff
md""" Se grafica un segundo de la muestra del audio Pink.ogg para mostrarla en detalle"""

# ╔═╡ 1e7678f0-73ce-11eb-02a6-053cc2e84b2f
plot(tseg,tst[n],
label = "f = 44100",
legend = true,
xlabel = "Tiempo(Seg)",
ylabel = "Amplitud",
#xlims=(0, 1500),
#ylims=(-1.1, 1.1)
)

# ╔═╡ 202f4c80-73ce-11eb-2b03-1dc6c7240993


# ╔═╡ b9ad22ac-3e67-11eb-35e1-7f4579b64838
md"""
#### Ejercicio 2)

**Mediante la función `fft`, realice un análisis espectral de la señal mostrando la densidad espectral de potencia (PSD) en función de la frecuencia en Hz.  Verifique que la mayor PSD se concentra en las bajas frecuencias. Utilice escala logarítmica para mostrar esultado.**

**La PSD, ``S_{XX}``, la puede estimar como ``S_{X X}(\omega) = 1/T |X(\omega)|^2``, donde ``T`` es la duración del segmento temporal que usa para la estimación, o mejor aún, dividiendo al audio en segmentos, realizando la estimación anterior con cada uno, y finalmente promediándolas.**
"""

# ╔═╡ 2b430db0-7963-11eb-3585-893f8168af1c
md"""Para realizar lo pedido se hace uso de la función **fft()** quien recibe las muestras y calcula la _transformada discreta de Fourier_, para calcular la PSD según lo indica la formula dada se necesita la continua por ello habrá que dividir por la frecuencia de muestreo.
Luego la función hará uso de ellos calculándola según lo sugerido, se divide el audio en segmentos obteniendo la PSD para cada uno y promediándolas
."""

# ╔═╡ 4fc8c804-3ea4-11eb-3e97-eb6709f1c0aa

## retorna el PSD en dB.

function PSD(x, n, t_win)
	
	return mean(10log10.(abs.( fft(tst_seg) ./ sr ).^2 ./t_win)
	for tst_seg in arraysplit(x, n, 0)
	)
end

# ╔═╡ 56f4c0b0-73ce-11eb-09fd-fba3b409e4d9
begin
	psd = let
			n = 1000; 			#4410 	#cantidad de muestras por segmento
			t_win = n / sr; 			#ancho del segmento.
			PSD(tst, n, t_win)
	end

	f_psd = let
		m = 1:length(psd)/2;
		(m.*sr)./length(psd)
	end

end

# ╔═╡ d987dc10-7621-11eb-26d9-d148f796fc11
md""" Se procede a graficar la PSD, notar que para las frecuencias más bajas es donde la PSD es mayor, esto sugiere que lo importante de la señal se concentra en un rango menor de frecuencias y se podría despreciar lo que sucede a más altas frecuencias. """

# ╔═╡ 7776221e-73ce-11eb-003c-857b22a7e7c1
plot(f_psd,
	psd[1:div(length(psd),2)],
	xlims = (0, 6000),
	ylims = (-80, -40),
	title = "PSD vs Frec.",
	xlabel = "Frecuencia (Hz)",
	ylabel = "PSD en dB (1/Hz)",
	legend = false
	)


# ╔═╡ bffcb902-73ce-11eb-31c9-05d7c512a3d0


# ╔═╡ b60ae59e-3e67-11eb-123e-11c0cba7d09e
md"""
#### Ejercicio 3)

**Mediante la función `fft`, obtenga el espectro de una ventana rectangular y una de Hamming de igual duración y grafique su potencia en escala logarítmica en un mismo gráfico para compararlos. ¿Cuál es la resolución en frecuencia de cada ventana? Al realizar un espectrogama, ¿qué ventaja tiene utilizar la ventana de Hamming frente a la rectangular?**
"""

# ╔═╡ de8b8c70-7964-11eb-3127-d3731697c053
md""" Se procede a graficar su densidad espectral de energía siendo que este aporta lo mismo para lo que se quiere analizar.
La curva azul corresponde a una ventana rectangular mientras que la roja a la ventana de hamming
 """

# ╔═╡ 4e904a84-3ea4-11eb-0c12-b1fccd5f7036
let
	normalize(x) = x ./ sum(x)
	plot()
	plotwin(win) = plot!(range(0; stop=2pi, length=5000 ),
		log.( ( abs.(fft(padright(normalize(win(41)), 5000#=501=#))) ).^2 );
		legend= false,
		xlims=(0, pi), ylims=(-15, 2), 
		ylabel=" DSP dB ",
		xlabel= " Ω ",
		title = "densidad espectral de energia")
	foreach(plotwin, [rect, hamming],)
	plot!()
end

# ╔═╡ 98aab180-73cf-11eb-092d-cb520eba3c36
md""" Se puede observar la diferencia entre los lóbulos secundarios entre la ventana rectangular y la de hamming, la relación entre el primer lóbulo y el segundo es mucho mejor con hamming.

Dado que el grado de fugas depende de la amplitud relativa del lóbulo principal frente a los lóbulos laterales, la ventana de hamming tiene mejor compromiso entre estas (a costa de un decaimiento más lento en alta frecuencia), además se sabe que la anchura del lóbulo principal y la amplitud relativa de los lóbulos secundarios dependen del ancho de la ventana.
 """

# ╔═╡ d61f59b0-73d1-11eb-0edb-938b1f3defa2
let
	order = 501;
	Hamm = specplot(tst; fs=sr, window=hamming(order), title= "Hamming")
	Rect = specplot(tst; fs=sr, window=rect(order), title= "rectangular")
	
	plot(Hamm, Rect)
end



# ╔═╡ f538d140-73d2-11eb-2248-abda76c139a5
#=md"""
El problema que tiene la ventana rectangular esta en sus discontinuidades, haciendo que lal señal se
distorcione, esto se puede ver mejor en la DSP donde claramente se ve como aparece un lobulo
secundario bastante proximo al principal a comparacion del de hamming donde el lobulo secundario
se vé mas atenuado dado que este no tiene discontinuidades. Si se observa los diagramas espectrales
se puede ver que en la ventana rectangular para frecuencias muy altas todavia hay picos de amplitud
elevados, mientras que en la de hamming estos "picos" se centran en las bajas freucencias"""
=#

# ╔═╡ c0621930-796e-11eb-1a94-a118feca6e6b
md""" La resolución en frecuencia tiene que ver con el ancho del lóbulo principal, esto porque al multiplicar la señal por una ventana se convolucionan en el espectro.


El ancho del lóbulo principal tiene un papel importante en la resolución en frecuencia. 
Si el ancho del lóbulo principal se estrecha es posible diferenciar dos componentes de frecuencia cercanas. Cuando el lóbulo principal se estrecha, con lo cual se incrementa la resolución en frecuencia, la energía de la ventana es distribuida entre los lóbulos laterales y las pérdidas espectrales son más grandes.

De manera ideal se busca que el espectro de una ventana deba parecer a un impulso y estar confinado a un lóbulo lateral tan estrecho como sea posible y con tan poca energía en los lóbulos laterales como sea posible

."""

# ╔═╡ 1ae45790-7a1b-11eb-03b9-a3212ee88756
md""" Según el _Oppenheim-Schafer_ se puede estimar el ancho del lóbulo principal de las siguientes formas :

- Hamming: $\frac{8 \pi}{M}$= ~0,62  		con M=41.
- Rectangular: $\frac{4 \pi}{M+1}$= 0,3  	con M=41.

Midiendo el ancho del lóbulo principal se obtiene 0.64 para Hamming y 0.28 para la rectangular siendo estas la resolución en frecuencia.

Notar que el ancho del lóbulo principal de la ventana de hamming es el doble que la rectangular.

"""

# ╔═╡ b2025250-3e67-11eb-39a2-73292bbf17c9
md"""
#### Ejercicio 4)
**Implemente un sistema para reducir la frecuencia de muestreo del audio, de 44100 Hz a 5512.5 Hz. Muestre un diagrama en bloques y justifique su diseño. Diseñe el filtro pasabajos mediante el método de ventaneo explicando y justificando las decisiones, de forma tal que su retardo sea menor a 1 ms. Graficar un diagrama de polos y ceros, respuesta en frecuencia (módulo y fase), respuesta al impulso, y retardo de grupo.**
"""

# ╔═╡ 4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
md"""
La velocidad de muestreo puede ser reducida por un factor entero M= 8, haciendo una declinación de la señal discreta generando una nueva secuencia $ x_d[n] $ .

$$x_d[n] = x_d[m*M] = x_c[m*MT]$$, con T: periodo de muestreo.

En el espectro esto se puede ver con amplitud 1/N y ancho de banda desde 
$-M*\Omega_m$ a $M*\Omega_m$
entonces la reducción en frecuencia de muestreo puede generar un traslape para $\Omega_N > \pi/M$ . 

Por ello para no perder información se realiza un filtro anti-aliasin antes del submuestreo.



"""

# ╔═╡ 3936fb00-73ed-11eb-25cd-812cdcf591fd
let
	img_path = "decimador.png"
	img = load(img_path)
end

# ╔═╡ c3a33cae-73d7-11eb-2707-d3604f680c0c
begin

fc= 2322#2756; 				#frecuencia de corte(hz) fc < fs/2 = (44100/8)/2.
OmegaC = (2*pi*fc)/sr;
end;

# ╔═╡ 1500ddb0-73e7-11eb-024b-1f68d5ac7525
md""" Observando el grafico de la PSD y teniendo en cuenta la cota dada por Nyquist de fc < sr2/2=~2756 se elige una frecuencia de corte razonable para la señal de 2322 para no estar tan al limite de Nyquist. """

# ╔═╡ b07a0860-79d3-11eb-197f-4ba2019e2521


# ╔═╡ 5fef73e0-73e7-11eb-36e6-d1a47f2c456c
md""" El retardo esta dado por (N - 1)/2 donde N es el ancho de la ventana, sabieno que las muestras y el tiempo estan relacionada por la frecuencia de muestreo, sera entonces:
 $\frac{N - 1}{2}$.$\frac{1}{sr}$ < 1mS. entonces para cumplir con lo pedido mi N < 88

"""

# ╔═╡ 27b2c6d0-73d8-11eb-33d5-e779bb4bf27d
begin
	H_ideal(Ω) = u(Ω + OmegaC) - u(Ω - OmegaC);
	
	plot(H_ideal, -π, π,
		xlabel= " Ω ",
		ylabel= "Amplitud",
	)
end


# ╔═╡ 79411da0-79ee-11eb-3f3d-5b00dcffd43c
md""" Se arma un filtro ideal y su respuesta al impulso."""

# ╔═╡ ece5f452-79ef-11eb-0b5f-8141791cbab8
md"""Grafico la respuesta al impulso del filtro ideal."""

# ╔═╡ 54934580-73d8-11eb-2211-3f5357c88358
#La respuestaal impulso de mi Hideal.
h_ideal(n, OmegaC) = (OmegaC/pi) * sinc(n *(OmegaC/pi)) #Omega = 2pi/fm

# ╔═╡ 340d0210-73d8-11eb-323c-e9501bdf682f
let ns = -40:40
	stem(ns, h_ideal.(ns,OmegaC))
end


# ╔═╡ fbde8120-79ef-11eb-05f7-2def33cc8263
md""" Grafico mi filtro convolucionado con la ventana con order 87."""

# ╔═╡ 7e90a16e-73d8-11eb-1e23-0b1e366467aa
#Grafico mi filtro ventaneado.
begin
	order = 87;
	ns = range(-(order - 1) / 2; length=order);
	winds = hamming(order);
	hs = h_ideal.(ns,OmegaC) .* winds;#<--filtro.
	stem(0:order-1, hs)
end


# ╔═╡ c4921f20-79ee-11eb-1cc5-99c64e998496
md"""Se crea el filtro mediante el metodo ventaneo con un order de 87 para tener un retardo menor a 1mS, ademas se la grafica en comparacion con la ideal"""

# ╔═╡ b18648f0-73d8-11eb-2f64-ff6fb23b096a
#Grafico el filtro obtenido con el pasabajos ideal.
begin
	samples = 500
	hfs = fft(padright(hs, samples))
	fs = range(0; step=2pi/samples, length=samples)
	plot(fs, abs.(hfs); xlims=(0, π))
	plot!(H_ideal, 0, π)
end

# ╔═╡ ff2934c0-79ee-11eb-05ea-b77dfd7f1a73
md""" De esta forma obtenemos un filtro FIR, estos filtros desplazan temporalmente a todas las frecuencia por un valor constante, es decir es de fase lineal, como se ve en el siguiente grafico."""

# ╔═╡ 92bed250-75e2-11eb-3e64-d15f6a219026
plot(fs, unwrap(angle.(hfs)); 
	xlims=(0, π), ylims= (-22, 0),
	xlabel= " Ω ",
	ylabel= "Fase",
	title= "Fase", label= false)

# ╔═╡ 9fee5cf0-75f8-11eb-0404-cd8ad839260c
md"""Dado que es fase lineal se procede a graficar el retardo de grupo del filtro, al ser la ase lineal la derivada se una constante.
El retardo de grupo se define como $-\frac{d(Fase(w))}{dw}$"""

# ╔═╡ 2fcff8a0-761d-11eb-12d0-1736a6b45e85
let
	fase= unwrap(angle.(hfs));
	fg = range(0; step=2pi/samples, length=samples-1)
	#retardo de grupo se define como -d(fase(w))/dw.
	plot(fg, -diff(fase)/(2π/length(fase)) ; 
		xlims=(0, π), ylims= (0, 50),
		xlabel= " Ω ",
		ylabel= "Retardo",
		title= "Retardo de grupo", label= false)
end

# ╔═╡ 29aa4a30-79f0-11eb-207d-6f7e985947d0
md""" Se procede a realizar un diagrama de polos y ceros """

# ╔═╡ 6ef5b922-73e3-11eb-0504-652233cd6dda
let 
	pr = polynomialratio(hs, [1]) 
	zplane(getzeros(pr), getpoles(pr))
	title!("Diagrama de polos y ceros")
end

# ╔═╡ d08961ee-7550-11eb-0529-513ae044a205
md""" Un sistema FIR es finita en el tiempo y esta definida por sus ceros, es decir no tiene polos en el plano finito salvo z= 0 (es decir para w->-inf).
Notar que al ser un pasabajos para Z= 1 se no hay ceros, ni en las cercanias, esto se corresponde a w= 0 y frecuencias bajas, es coherente con un pasabajos"""

# ╔═╡ 72acf4c0-73e3-11eb-1975-77a7c0a9c562
begin
m = 1000; #cantidad de muestras por segmento
t_win = m / sr;
#h[nsss]=2*fc*sinc(s*fc*nsss)
filt_tst = PSD(conv(hs, tst), m, t_win)
end


# ╔═╡ c3f1fca0-79f0-11eb-287e-4909d339ae38
md""" Se procede a graficar la PSD antes y despues de filtrar la señal con el filtro diseñado"""

# ╔═╡ 7c44b4a0-73e3-11eb-03c4-afb2354639a0
let
	#responsetype = Lowpass(fc; fs=44100)
	#designmethod = FIRWindow(hamming(500))
	plot(f_psd,
	psd[1:div(length(psd),2)],
	#xlims = (0, 6000),
	ylims = (-90, -40),
	title = "PSD antes y despues de filtrar.",
	xlabel = "Frecuencia (Hz)",
	ylabel = "PSD en dB (1/Hz)",
	legend = false
	)
	plot!(f_psd, filt_tst[1:div(length(filt_tst),2)])
end


# ╔═╡ cc11a150-73e3-11eb-0459-4fedbabe462d
md""" Notese que la diferencia que hay en los graficos es en escala logaritmmica por lo que luego de la frecuencia de corte la señal se atenua bastante, pero no es discontinia."""

# ╔═╡ 16be5cc0-7777-11eb-1b8e-b765ed6d37cd
begin
	
	function submuestrear(x,M)
		sample = x[1:M:end]
		return sample
	end	
	
#finalmente funcion que me submuestrea y devuelve e. .ogg a frecuencia de muestreo sr/8.
	function reduce_sampleRate(x,fm)
		M = 8 #44100 Hz a 5512.5 Hz
		return submuestrear(conv(hs,x),M)
	end
	
end

# ╔═╡ af4f3da4-3e67-11eb-3cc6-3378e0c12667
md"""
#### Ejercicio 5)

**Realice un espectrograma de la señal original y la filtrada y verifique los efectos del filtrado. Indique el tipo y longitud de ventana utilizada.**
"""

# ╔═╡ 3aa5434e-3ea4-11eb-20aa-b15564d4eb90
let
	order= 1024
	snf = specplot(tst; fs=sr, window=hanning(order), title= "sin filtrar" )
	sf  = specplot(conv(hs, tst); fs=sr, window=hanning(order), title= "Filtrada" )
	#sfg  = specplot(conv(hs, tst); fs=sr, window=hanning(order), title= "Filtrada" )
	plot(snf,sf, legend= false)
end

# ╔═╡ 38c05a80-73e4-11eb-274d-d538d2e3fb65
md""" Se observa que las muestras filtradas efectivamente tienen sus componentes de mayor amplitud por debajo de la frecuencia de corte luego para frecuencias mayores disminuye bastante pero no de forma discontinua."""

# ╔═╡ 982538c4-3e67-11eb-229e-dd2531a540d6
md"""
### Espectrograma

Los parámetros con los que se realice el espectrograma tienen una influencia directa sobre la extracción de las características. La elección de la ventana (tipo y longitud) es el parámetro principal, ya que debería lograr una resolución razonable tanto en tiempo como en frecuencia que permita distinguir las características espectrales y temporales con claridad.

#### Ejercicio 6)

**Realice un espectrograma de la señal sub-muestreada y muestre una imagen (con zoom) de alguna región donde se vean las características espectrales de la señal dada la ventana utilizada. Muestre y justifique cómo cambia la visualización de las características con 3 diferentes longitudes de ventana y comente qué longitud utilizará en el algoritmo. Realice las comparaciones con ventana rectangular y de Hamming.**
"""

# ╔═╡ f3971b40-73ea-11eb-3f2e-dd3a4443b9ab
ssm = reduce_sampleRate(tst, sr);  
#length(tst) #441000
#length(ssm) #55154 muestras
#441000/55154 =7.995793596112702->8

# ╔═╡ f906add0-74b1-11eb-156d-9b7c54bfd2ba
length(ssm)

# ╔═╡ 39f5fc86-3ea4-11eb-37f3-25feb7d2aee6
begin	
	order1= 2048  #orden de mis filtros
	order2= div(order1, 2)
	order3= order1*2
	
	sr2 = sr/8;   #Nueva frecuencia de muestreo.
	ovrlp=0.5#31/32;  #overlap propuesto por el paper
	#bajarle el overlap paraverlo mas suave	!!!!!	
end

# ╔═╡ f7de1ba0-79fc-11eb-1fa8-c9f36f88d7ba
md""" Se proceden a graficar los espectrogramas para cada el orden 2048 es se deduce experimentalmente viendo los gráficos.
también podemos estimarla con ayuda del paper dado para el TP, en donde proponen un ancho de ventana de 0,37seg. y dado que nuestra frecuencia de muestreo reducida es ~5512Hz se obtiene un ancho en muestras de 2040 muestras, ajustando un poco llegamos a las 2048
"""
#=Lafuncion plot no es muy amiga del specplot, se rompe a cadarato, por ello cargar imagenes de los espectrogramas y subirlas.
let
	#t_s=1/2322;
	t_win2=0.37 						#segun elpaper
	fm=5512.5
	muestras= t_win2 * fm 				#2040muestras.
		
end=#

# ╔═╡ e1906f20-7a00-11eb-3204-2f8cac355516
md""" el siguiente grafico muestra el espectrograma con el orden 2048 el cual se usara para el resto del trabajo."""

# ╔═╡ 5cc77020-74bc-11eb-154a-277125d7a831
let
	
	plot(
		specplot(ssm;fs=sr2,overlap= 31/32#=ovrlp=#,window=hamming(order1),title= "hamming"),
		specplot(ssm;fs=sr2,overlap= 31/32#=ovrlp=#,window=rect(order1),title= "rectangular"),
		layout= (1, 2),
		xlims= (1, 5),
		ylims= (0, 1500)
	)
end

# ╔═╡ fd9d3860-7a00-11eb-2d95-b5700362c015
md""" Si se grafica el mismo espectrograma pero con un orden 2 veces mayor, si bien a mayor orden mayor resolución espectral para este trabajo no es necesario dado que impondría mucha carga al programa y gráficamente se ve que tampoco es mucho mejor"""

# ╔═╡ aa6bbab2-7620-11eb-36d0-9b82996620bb
let
	order2= 1020
	order3= order1*2
	p1= specplot(ssm;fs=sr2,overlap= ovrlp,window=hamming(order3),title= "hamming");
	p2= specplot(ssm;fs=sr2,overlap= ovrlp,window=ones(order3),title= "rectangular");
	
	plot(p1,p2, layout= (1, 2),
		xlims= (1, 5),
		ylims= (0, 1500)
	)
end

# ╔═╡ 987a7410-7a01-11eb-1f48-dd5f50a1cd81
md""" Se grafica los espectrogramas con la mitad del orden utilizado, notar que disminuye la resolucion espectral."""

# ╔═╡ b52628a0-7620-11eb-1dda-65bc93222c5c
let
	order2= 1020
	order3= order1*2
	p1= specplot(ssm;fs=sr2,overlap= ovrlp,window=hamming(order2),title= "hamming");
	p2= specplot(ssm;fs=sr2,overlap= ovrlp,window=ones(order2),title= "rectangular");
	
	plot(p1,p2, layout= (1, 2),
		xlims= (1, 5),
		ylims= (0, 1500)
	)
end

# ╔═╡ 6d1698d0-74a7-11eb-15bf-35d2eae0061f


# ╔═╡ 9309e284-3e67-11eb-1ab2-612f6c748c3b
md"""
### Extracción de características

La función provista `stft` permite obtener la matriz $s$ con las sucesivas DFT de los segmentos de la señal ventaneados -- lo que se grafica en un espectrograma.

Estas DFT nos devuelven muestras del espectro en frecuencias equiespaciadas en escala lineal – junto con los vectores de tiempos y frecuencias correspondientes a cada columna y fila. 

A partir de esto, necesitaremos la energía en bandas de frecuencia equiespaciadas en escala logarítmica (esto es similar a la escala psicoacústica Bark) debido a que así funciona el oído humano, que es un buen reconocedor de canciones – la relación entre frecuencias y posición de resonancia de la membrana basilar es aproximadamente logarítmica.

Debemos obtener una matriz de energías `E`, cuyas columnas, al igual que las de la matriz del espectrograma `s`, representan características espectrales del $n$-ésimo segmento de señal ventaneado (frame `n`). Para todo `n`, el elemento `E[m, n]` de nuestra nueva matriz debe contener la energía total de todos los coeficientes de `s[k, n]` asociados a frecuencias que caen dentro de la $k$-ésima banda de frecuencia.

#### Ejercicio 7)

**Divida al espectrograma en 21 bandas de frecuencias, equiespaciadas en escala logarítmica (es decir, el cociente entre frecuencias de bandas consecutivas debe ser constante). La frecuencia inferior de la primera banda debe ser 300 Hz y la superior de la última, 2 kHz.** 
"""

# ╔═╡ e7793b50-7a01-11eb-35ce-bfc0ff756a76
md""" Se define un vector separado logarítmicamente entre 300 y 2k, con el que se trabaja en adelante para dividir el espectrograma en 21 bandas."""

# ╔═╡ 5f636b02-3ea4-11eb-3f78-6f693a936992
fbands= exp.(range(log(300); stop=log(2e3), length=22))

# ╔═╡ 59fddd50-7513-11eb-1f17-07ef0c9f45f0
md""" cada frecuencia es igual al anterior multiplicada por el mismo factor """

# ╔═╡ 6813ab40-7513-11eb-11b4-8f3040b56c91
fbands[2: end] ./ fbands[1: end - 1] 		#coeficientes constantes.

# ╔═╡ 22d04590-7a02-11eb-3fe2-ab6b63005994
md""" Esra banda va de 300Hz a 2kHz, tiene 22 frecuencias para 21 bandas"""

# ╔═╡ a46234e0-7513-11eb-21d2-8592c8fec31d
fbands[[1, end]]  				#frecuenciasinferiores  y superiores acordes.

# ╔═╡ aeb41d00-7513-11eb-0980-9b6250c020af
length(fbands)

# ╔═╡ 61106600-7a02-11eb-256f-d724de274557
md"""Se genera una matriz **S** quien posee las sucesivas DFT dada por la función stft() en donde se tomaran la primer mitad de las filas."""

# ╔═╡ e6a127a0-751b-11eb-34ac-418d1ad01639
begin
	ovlp= 1827
	#floor(Int, (31/32)*length(hamming(order1)))
	#2048*(31/32) fui bajando el overlap para que me quede mejor la huella.
	s= stft(ssm; overlap= ovlp, window= hamming(order1), nfft= order1);
	S=s[1:div(2048, 2), :];
end

# ╔═╡ d3fe3e30-7a02-11eb-25ea-85d66902a047
md""" Luego entonces se genera un vector similar con las posiciones correspondiente a cada frecuencia del vector de frecuencias definido anteriormente, los cuales se necesitan para operar con la matriz **S**"""

# ╔═╡ 46b25e4e-7609-11eb-3ef8-bff9a5a5898d
begin
	#fil, col = size(s);
	fil, col = size(S)
	div_fil_f = floor.(Int, (fil.* fbands[1:end]) ./ (sr2/2))
end


# ╔═╡ 7fda0ee0-752c-11eb-2b7d-5bba9b2a3592
#=
begin
	function mean_nfil(mat, i, j, wfil)
		@assert wfil >= 1
		aux=(abs(mat[i,j]))^2;
		[aux = ((abs(mat[x,j]))^2 + aux) for x=i+1:1:i + wfil -1]
#aca quien carajos es x!!!
		return aux
	end
#recibe la matriz "mtr" su fila inicial y el ancho de filas quedebe promediar.

	function mean_band_gap(mtr, init_fil, fil_band)
		fil_mtr, col_mtr = size(mtr);
		mat2= mtr[1, :];

		mat2= [mean_nfil(mtr, init_fil, j, fil_band) for j= init_fil:1:col_mtr]

		return mat2
	end
end

=#

# ╔═╡ 8a696ed0-7524-11eb-0764-214ebea3d1e7
#=
E = let
	
	fil_band(w_seg) = sr2*(fbands[w_seg] - fbands[w_seg + 1]); ##muestras que hay en 	      															esas frec
	fil_mtr, col_mtr = size(s);
	XD = zeros(length(div_fil), col_mtr)
	
	XD= hcat([ mean_band_gap(s, i , w_seg)
				for w_seg in div_fil[2:1:end] 
				for i in div_fil[1:1:end-1]
				])
				#for filE in 1:1:length(div_fil) ]
	
end
#imposible!!!! no puedo no puedo!!! quiero armar con los vectores una matriz!!!!
=#

# ╔═╡ 0fa90dc0-7a03-11eb-3301-5166ff8dfc61
md""" Finalmente se genera una función que para cada fila "m" y columna "n" retorna el valor de la energía correspondiente a una matriz de energía generada a partir de la matriz **S** """

# ╔═╡ 4a342420-75df-11eb-35be-0bc9cc550163
#para cada banda m en s calculla la media de los abs ^2 de las filas correspondientes.
E(m,n)= mean( abs.( S[div_fil_f[m]:div_fil_f[m+1], n] ).^2 )

# ╔═╡ 8deaf928-3e67-11eb-0327-31e0f74de814
md"""
Finalmente, las características se obtienen mediante una función que opera sobre E(m, n) según:


![sys-tpe4.png](https://i.imgur.com/dGoWXWE.png)

En palabras, $H[m, n]$ es 1 si la diferencia de energía entre las bandas $m$ y $m+1$ para el *frame* actual n es mayor a la del *frame* anterior $n-1$. Experimentalmente, se verficó que estas características son robustas ante varios tipos de procesamientos y distorsiones del audio\[5\].

#### Ejercicio 8)

**Implemente el algoritmo de extracción de características, calculando la huella. Debido al efecto borde, $H[m, n]$ debería resultar una matriz de 20 filas. Ejecute el algoritmo sobre un segmento de audio y muestre una imagen de la huella digital acústica obtenida.**
"""

# ╔═╡ 198f8280-7a06-11eb-2a86-d16bc498aaef
md""" Para construir la **H**, la matriz huella, se hace uso de la función _E()_, generando una subfunción _bern()_ que según el criterio de diferencias de energías mostrado en la definición de **H[m,n]** retornara _1_ si la condición se cumple _2_ en caso contrario.
Luego entonces a partir de la subfunción mencionada se genera una función _H()_ que tornara la matriz de huella.

"""

# ╔═╡ 6476e9fc-3ea4-11eb-3873-b765108f4bab
begin
	
	function bern(m, n)
	  (E(m+1, n)-E(m,n)) > (E(m+1, n-1)-E(m, n-1)) && return 1
	  return 0		
	end 

	function H(s)
		fils, cols = size(s)
		filE= length(div_fil_f)-2
	
		#[H[m,n] = bern(m, n+1) for m in 1:1:filE for n in 2:1:cols-1]
		H = [bern(m, n+1) for m in 1:1: filE,  n in 1:1:cols-1]
		
		return H
	end
end

# ╔═╡ 2435f260-7618-11eb-1b22-d3eda517d5f0
md""" 
Usando la función H() se genera a modo de prueba la huella para las muestras de "Pink.ogg"
"""

# ╔═╡ d33d59f0-75f6-11eb-2c74-1b8e3be2f9c1
huella= H(S)

# ╔═╡ feb5ec9e-7a06-11eb-034a-c169e0377755
md""" 
Para verlo mejor se grafica la matriz en escala de grises pudiendo ver una imagen con la huella obtenida.
"""

# ╔═╡ 5bf232a0-7568-11eb-302c-a17935859a8a
plot_huella(h)= Gray.(float.(h));

# ╔═╡ 42c730b0-7625-11eb-1c4e-6b516c1eb2e7
plot_huella(huella)

# ╔═╡ 7b0ea820-7771-11eb-2c41-15bf560369ab
#=
begin
	
	sp= stft(x5k, window= hamming(2048), nfft= 2048, overlap= 1827, fs=sr5k)
	sppow= abs.(sp).^2
	freqs= range(0:step=sr5k/size(srpow, 1), length= size(srpow, 1));
	
	nbans= map(fr ->findfirst(fr .< fbans), freqs)
	nrg= zeros(length(fbans)-1, size(srpow,2))
	
	for i in 1:length(fbands)-1
		nrg[[i], :] .= sum(sppow[nbans . == i+1, :]; dims=1)
	end
	
	huella= diff(diff(nrg; dims=2);dims=1) .>0
end
=#

# ╔═╡ 89743a62-3e67-11eb-209e-9b1f3cc84e34
md"""

## Confección de la base de datos

En principio, la base de datos simplemente debe

* Guardar, para cada *frame*
  * las características obtenidas
  * un identificador de la canción de la cual provino
  * el número de *frame* dentro de la canción.
* Permitir buscar el identificador de canción y el número de *frame* a partir de una característica, o informar si la característica buscada no se encuentra en la base de datos.

Así, cuando se desea identificar una música desconocida, se calculan las características de cada *frame*, se obtienen los identificadores de la base, y se opera con ellos de algún modo razonable para devolver la (o las) canciones más probables.

Ahora bien, más allá de este trabajo, es importante que este tipo de algoritmos escalen a bases de datos grandes y funcione rápido; para eso, es crucial que la la base de datos aproveche bien el espacio en memoria y permita una búsqueda eficiente. Por esto, se le proveen funciones que implementan una versión sencilla de una base de datos que cumple con estas características, que deberán ser entendidos.

Para confeccionar la base de datos, se utilizará la función *generar_DB*. Esta función requiere que el alumno haya definido otra función que, dado un archivo de audio como entrada, lo pre-procese, analice y extraiga su huella digital acústica en forma automática. 


#### Ejercicio 9)

**Encapsule su algoritmo de extracción de huellas acústicas en una función llamada `generar_DB` que reciba como entrada el path del archivo de audio y devuelva su huella digital acústica.**

```julia
\"""
	generar_huella(fname)

Devuelve la huella del audio en el archivo `fname`.

Tiene que estar muestreado a 44100 Hz.
\"""
function generar_huella(fname::String)
	#...

	return huella
end
```
"""

# ╔═╡ d4bbda22-7797-11eb-1835-ddf37a08cc0f
md"""
Como indica el enunciado se procede a juntar todas las funciones definidas anteriormente en una para generar la huella y así poder usarla función en cualquier otra canción en muestreada a 44100Hz.
Para ello se redefine la función _E_ y _bern()_ creando las funciones _EE_ y _bern2_ las cuales permiten ser usadas más elegante y eficientemente en la función generar_huella().

"""

# ╔═╡ 304f1ef0-7776-11eb-1926-13a4a3e50e26


# ╔═╡ 76604f60-777e-11eb-2543-09be49dcc4b1
#Gnero las mismas funciones de antes modificadas para poner usarlas mejor en generar_huella.

begin
#calcula la energia para cada fila y columna de la matriz s.	
	function EE(s, nb, m, n)
		return mean( abs.( s[ nb[m]:nb[m + 1] , n] ).^2 )
	end
#devuelve 0 o 1	segun se cumpla la condicion.
	function bern2(s, nb, m, n)
	  (EE(s,nb,m+1,n) - EE(s,nb,m,n))>(EE(s,nb,m+1,n-1)-EE(s,nb, m, n-1)) && return 1
	  return 0
	end 
end

# ╔═╡ bd363c20-7a07-11eb-1a07-a1e3ea91852f
md"""
La función _generar_huella()_ debe recibir una canción, pero más adelante esto traería problemas al tener que convertir muestras a canción y luego canción a muestras, por ello se define otra subfunción _generar_huella!(), que hace lo mismo que la anterior pero recibe las muestras de la canción a procesar.
Entonces generar_huella() usa a generar_huella!() convirtiendo la canción en muestras previamente.

"""

# ╔═╡ 4b0b7360-7776-11eb-1bd5-ddd0c82e9e06

#genera la huella de x.	
function generar_huella!(x::Vector) 				    	#debe recibir un vector.
	
	xsrsub = reduce_sampleRate(x, sr); 					  #sr=44100 variable global.
	rsr = sr / 8

	order= 2048
	ovlp= 1942 		#Segun lo calculado en el ejercicio 10)           #1827

	s= stft(xsrsub; overlap= ovlp, window= hamming(order), nfft= order);
	S= s[1:div(order, 2), :] 			#se puede hacer en una linea??

	fbands= exp.(range(log(300); stop=log(2e3), length=22))   #divido en bandas.
	filS, colS = size(S);
	div_fil = floor.(Int, (filS .* fbands[1:end]) ./ (rsr/2)) #busco el indice.

	filE= length(div_fil)-2

#devuelve la matriz H.	
	H= [bern2(S, div_fil, m, n + 1) for m in 1:1: filE,  n in 1:1:colS-1]

	return H
end;
	


# ╔═╡ d4179090-77a7-11eb-3c1b-4da6e8152c16
function generar_huella(fname::String )	
	generar_huella!(to_mono(fname))
end

# ╔═╡ 5795e360-7a08-11eb-0142-31c7502d8428
md"""
A modo de prueba se genera nuevamente la huella para _Pink.ogg_ viendo que es similar a la anteriormente mostrada solo que esta vez fue generada a partir de la función generar_huella().
"""

# ╔═╡ 741304fe-3ea4-11eb-15e8-09908d98ecb3
hhuella= generar_huella("Pink.ogg");

# ╔═╡ ed9dce70-7785-11eb-27b2-8d22fa913c96
plot_huella(hhuella)

# ╔═╡ 855a7d2e-3e67-11eb-0f46-a5c786d5caf3
md"""

#### Ejercicio 10)

**Observe que la cantidad de elementos a guardar en la base de datos se incrementa conforme la longitud de las ventanas del espectrograma inicial disminuye, o el solapamiento entre ventanas se incrementa. Determine el solapamiento entre ventanas del espectrograma para obtener una densidad de aproximadamente 25 elementos por segundo y utilice este valor para el ejercicio siguiente.**
"""

# ╔═╡ 4df7d760-7796-11eb-2563-932eb8c9fece
md""" 
Mediante la siguiente función se calcula la densidad pedida, de manera de encontrar 25 muestras por segundo según las columnas de la matriz h, luego con ello se encontró que con un overlap de 1942 y una ventana de 2048 se consigue lo pedido.
"""

# ╔═╡ 7a02bdf0-7793-11eb-3606-9168a11719a4
#col(H)/(tiempo de la canción[seg])
function density(length_song, colH, fr)
	
	time_song= length_song / fr
	dst= colH / time_song
	return dst
end

# ╔═╡ c19daa60-7795-11eb-03a6-d134602f316a
filH, colH = size(hhuella)

# ╔═╡ e66049e2-7793-11eb-0560-1189121d97f9
density(length(tst), colH,sr/2)  #muestras por segundo

# ╔═╡ 81717fc8-3e67-11eb-05fc-5bde46597f8a
md"""

#### Ejercicio 11)

**Ejecute la función `generar_DB` para confeccionar la base de datos completa de su lista de canciones. Utilice al menos 40 canciones para llenar la base de datos. Puede usar la lista de canciones provista, y/o usar una lista de canciones propia. (Recuerde verificar que la frecuencia de muestreo de sus canciones sea de 44100 Hz).**
"""

# ╔═╡ 8a2d6fe0-7a09-11eb-1946-694ee29e9d86
md""" 
Se cargan las muestras del subdirectorio _40songs_ y se usa la función _generar_DB_ provista para generar la base de datos, el cual se utilizara para comparar una dada huella  con las que se encuentran en la base de datos y así poder encontrar la información de la canción que el usuario este buscando."""

# ╔═╡ b91537ac-3ea4-11eb-14d6-d341c535d83e
begin
		# Cambiar por el nombre del subdirectorio donde estén sus canciones
	songsdir = "40songs"	
	songs = readdir(songsdir)
end

# ╔═╡ 73333e92-3e85-11eb-26b6-7f0309ef2ee9
@bind songpicked Select(songs)

# ╔═╡ feb5d512-3e85-11eb-0116-29e4d9539595
LocalResource(joinpath(songsdir, songpicked))

# ╔═╡ 0cf7ba9c-3e74-11eb-18e2-c38aa20f9e9a
# Código de `generar_DB`
begin
	
	"Guarda la información para identificar un frame"
	struct FrameID
		songID     :: UInt16
		frameindex :: UInt32 # índice dentro de la base de datos
	end
	
	"Genera la base de datos a partir de la lista de archivos de audio"
	function generar_DB(songs; dir="")
		
			# Inicializa la base, que es un vector de 2^20 vectores de FrameIDs
		db = [ copy(FrameID[]) for _ in 1:2^20]

			
		for songid in 1:length(songs) # Para cada canción
			
				# Genera la huella
			huella = generar_huella(joinpath(dir, songs[songid]))
			
				# La agrega a la base de datos
			for i in 1:size(huella, 2) # Para cada frame en la huella de la canción
					# Agregar el landmark (huella del frame) a la base
				addlandmark!(db, huella[:, i], songid, i)
			end
		end

		return db
	end
	
	
		
		# Agrega el landmark (huella de un frame) a la base de datos, asociado a la
		# canción con ID `songid` y número de frame `frameidx`.
	function addlandmark!(db, lm, songid, frameidx)
		fr = FrameID(songid, frameidx)
		idx = landmark_to_index(lm)

		push!(db[idx], fr)
	end
	
		# convertir un landmark -- vector de 20 bools -- en un entero 
		# con esa representación binaria
	landmark_to_index(h) = sum(2 .^ (0:19) .* h) + 1
end;

# ╔═╡ 928af020-7a0b-11eb-325c-81c84be31414


# ╔═╡ e9255b8c-3e74-11eb-2960-5d01b0c99b13
db = generar_DB(songs; dir=songsdir);

# ╔═╡ 5d0a5800-779a-11eb-2f2c-335d6f94ef06
db
#db es un vector de vectores!

# ╔═╡ 7c7c1424-3e67-11eb-1da0-5dbad0171b20
md"""
# Test del algoritmo

En esta etapa final, pondremos a prueba la eficacia del algoritmo completo de reconocimiento. El procedimiento consistirá en obtener el porcentaje de aciertos del método al reconocer segmentos de audio, los cuales serán sometidos a distintas distorsiones.

Para esto se suministra la función *query_DB*, la cual recibe como entradas la base de datos y la huella acústica del segmento de audio a reconocer, y devuelve el ID del la canción que mejor coincide con la huella. Para entender cómo opera esta función, lea el Apéndice y el código.
"""

# ╔═╡ 415e32e6-3e76-11eb-17fa-23bd653fb975
function query_DB(db, huella::AbstractMatrix)
	
		# Cantidad de landmarks en la huella
    n = size(huella, 2)
		
		# Extraigo la información de todos los frames que encajaron con algún landmark
		# Vea que puede haber más frames que landmarks si hubo firmas repetidas
    frames = [ fr for lm in eachcol(huella) for fr in db[landmark_to_index(lm)]]

		# Función que toma un songID y devuelve un puntaje según cuántas veces aparece
		# como máximo en `frames` en algún intervalo de tiempo de `n`
    function score(sid)
			
			# Rescato los índices de frame asociados al songID pedido
        fids = [ fr.frameindex for fr in frames if fr.songID == sid ]
		
			# Armo un vector con deltas de peso 1 ubicadas en los índices de frames
        aux = zeros(Int, maximum(fids))
        aux[fids] .= 1
		
			# Convoluciono eso con una ventana de ancho `n` y devuelvo el máximo
        return maximum(conv(ones(Int, n), aux))
    end

		# Listo todos los songIDs candidatos, borrando duplicados
    sids = union(getfield.(frames, :songID))

		# Devuelvo el songID que maximice el `score`
    return sids[findmax(score.(sids))[2]]
end;

# ╔═╡ 76ce23dc-3e67-11eb-0be0-91b6781840fb
md"""

#### Ejercicio 12)

**Evalúe la tasa de aciertos del algoritmo identificando segmentos de duración $T$ con tiempo inicial elegido al azar de canciones elegidas al azar (vea la función `rand`). Las canciones deberán ser las mismas que utilizó para confeccionar la base de datos. Realice la evaluación para 50 segmentos con duración T  entre 5, 10 y 20 segundos cada vez (150 evaluaciones en total) obteniendo la tasa de aciertos en cada caso.**
"""

# ╔═╡ a9b2ac20-7a0b-11eb-13cd-2de6007fc26d
md""" 
Para evaluar la tasa de aciertos se usa la función _query_DB()_ en donde para un segmento de "t", tiempo dado se selecciona un canción de forma aleatoria y se le obtienen las muestras a partir de un tiempo inicial aleatorio, luego se le genera la huella a dicho segmento y con _query_DB()_ se hace la búsqueda en la base de datos.
si este segmento de canción se encuentra en la base de datos retorna _true_ en caso contrario _false_
.
"""

# ╔═╡ c970df70-77ac-11eb-384a-9b8a5e2e9547
md""" 
Para calcular la tasa de acierto en los tiempos indicado se hace el promedio entre 50 muestras para cada tiempo y será retornado como vector el cual contendrá la media de aciertos correspondientes."""

# ╔═╡ 00da49d0-7602-11eb-31f1-dfe5b07d6349


# ╔═╡ 7229577a-3e67-11eb-0c71-f383056175d1
md"""
#### Ejercicio 13)

**Repita el ejercicio 12 sumando ruido a los segmentos de audio. Utilice la función `randn` para generar las muestras de ruido. Evalúe tasa de aciertos para $SNR =0dB$, $10dB$ y $20dB$, mostrando sus resultados en una tabla para 9 combinaciones de longitud temporal y ruido. Nota: $SNR=10 log_{10}(P_X/P_N)$ donde $P_X$ es la potencia media de la señal sin ruido, y $P_N$ es la potencia media del ruido sumado a la señal. Para el cálculo de la potencia media puede utilizar la función `var`, que estima la varianza de una señal, ya que las señales de audio no deberían componente continua o valor medio.**
"""

# ╔═╡ 746ccbb0-7a0e-11eb-11f0-f107a2b47442
md""" 
Para hacer las pruebas con ruido agregado se modificó la función genera en el ejercicio anterior agregando en los argumentos el parámetro "SNR" el cual es cero por defecto, y un parámetro booleano el cual me activa un segmento de Código que agrega el ruido, "noise= false" por defecto, poner en _true_ para agregar ruido.
"""

# ╔═╡ 841f7430-77b2-11eb-180b-11a3ec034424
function add_noise(x, snr)
	xn= x 									#copio la señal original.
	snr= 10^(snr/10)  						#snr_db=10log10(snr)
	#potencia media es suma de los valores ^2 /largo
	pn = var(xn)/snr       				 	#el ruido y la señal tienen media nula.
	
	#xn= [x .+sqrt(pepe) * randn(length(x)) for i in 1:length(x)]
	xn= x .+sqrt(pn) * randn(length(x))    #genero y sumo ruido
	
	return xn
end

# ╔═╡ 7d82ea60-779b-11eb-2e75-8997f2c88b53
#toma una duracion aleatoria procesa la huella y devuelve true si fue encontrado exitosamente.
function is_song(t; snr=0, noise::Bool=false)#keyword arguments
	
	nduration= floor(Int, t * sr )					#seg-->muestras.
	song= rand(songs)								#elijo una cancion al azar.
	rsong= to_mono(joinpath(songsdir, song))		#cargo la cansion.
	ninit= rand(1: length(rsong) - nduration) 		#elijo inicio al azar.
	x= rsong[ninit:ninit + nduration - 2] 		    #obtengo el segmento a testear.
	
	if noise == true
		x= add_noise(x, snr)
	end
		
	
	return song == songs[query_DB(db, generar_huella!(x) )]
end

# ╔═╡ b99bea60-7601-11eb-3f06-ef522296b39d
 [mean(is_song(t) for _ in 1:50 )
		for t 	in (5, 10, 20)]

# ╔═╡ 8fe96ca0-77ce-11eb-1d58-27201d27fc2f
md""" 
Luego se calculan para 9 tiempos predefinidos las  tasas de aciertos en las mismas condiciones del ejercicio anterior.
"""

# ╔═╡ 129267b0-77b6-11eb-2fa3-5914fd885d43
#esto se toma su tiempo, descomentarlo para las pruebas, mejor tiempo 22minutos 50 segmentos.
[mean(is_song(t; snr= algo, noise=true) for _ in 1:50 )
		for t in (3.5, 7, 12, 14.9, 17, 15.6, 23.5, 24, 24.7), 
			algo in (0, 10, 20)
		]

# ╔═╡ 7b6294de-77bb-11eb-2111-01df80d418e5
#[filas*columnas for filas in (5, 10, 20), columnas in (0, 1, 2)]

# ╔═╡ 713783c0-77bd-11eb-2b34-17a3651f41c9
md""" 

a continuacion se muestras distintas tablas las cuales contienen las tasas de acierto  con ruido agregado de 0dB, 10dB y 20dB.

Tabla para 50 segmentos con ruidos para t= 5,10 y 20 segundos.

|T[seg]| 0dB    | 10dB   | 20dB |
|:-:   | :----: | :----: | ---: |
|5     | 0.94   | 1.0    | 1.0  |
|10    | 0.98   | 1.0    | 1.0  |
|20    | 0.98   | 1.0    | 1.0  |



Matriz para 50 segmento con ruidos, para t= 4.5, 7, 12, 15, 17, 20, 23, 25 y 27 segundos.


|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5 |0.86  |1.0  |1.0
|7   |0.96  |1.0  |1.0
|12  |0.98  |1.0  |1.0
|15  |0.98  |1.0  |1.0
|17  |1.0   |1.0  |1.0
|20  |1.0   |1.0  |1.0
|23  |0.98  |1.0  |1.0
|25  |1.0   |1.0  |1.0
|27  |1.0   |1.0  |1.0

Matriz con ruidos para 50 segmentos de t: 3.5, 7, 12, 14.9, 17, 15.6, 23.5, 24y 24.7 segundos.

|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|3,5 |0.88  |1.0  |1.0
|7   |0.92  |1.0  |1.0
|12  |1.0   |1.0  |1.0
|14,9|0.98  |1.0  |1.0
|17  |0.98  |1.0  |1.0
|15,6|1.0   |1.0  |1.0
|23,5|1.0   |1.0  |1.0
|24  |1.0   |1.0  |1.0
|24,7|1.0   |1.0  |1.0

|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|3,5   |0.92 |1.0 |1.0
|7     |0.98 |1.0 |1.0
|12    |0.98 |1.0 |1.0
|14,9  |1.0  |1.0 |1.0
|17    |0.98 |1.0 |1.0
|15,6  |1.0  |1.0 |1.0
|23,5  |0.98 |1.0 |1.0
|24    |1.0  |1.0 |1.0
|24,7  |1.0  |1.0 |1.0

Notar que en todas las mediciones se obtuvieron aciertos de ~100% ecepto para el ruido.
Dado que con 0dB es cuando hay mas ruido sera mas complicado generar una huella que pueda ser encontrada en la base de datos, aun asi la tasa de aciertos para 0dB suele ser buena, mayor al 80%.


"""

# ╔═╡ 612586c0-7a25-11eb-247f-8f55aa84cb91
add_noise(tst,0)|>sound

# ╔═╡ 8274ebc0-7a13-11eb-2632-3f936f3eb7e2
#=
|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5|1.0|  1.0|  1.0|
|7  |1.0|  1.0|  1.0|
|12 |1.0|  1.0|  1.0|
|15 |0.0|  1.0|  1.0|
|17 |1.0|  1.0|  1.0|
|20 |1.0|  1.0|  1.0|
|23 |1.0|  1.0|  1.0|
|25 |1.0|  1.0|  1.0|
|27 |1.0|  1.0|  1.0|

La siguiente tabla muestra las tasas de acierto para 5 segmentos con duracion de segmentos t = 4.5, 7, 12, 15, 17, 20, 23, 25 y 27 correspondientes.

|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5|1.0  |1.0  |1.0
|7  |1.0  |1.0  |1.0
|12 |1.0  |1.0  |1.0
|15 |1.0  |1.0  |1.0
|17 |1.0  |1.0  |1.0
|20 |1.0  |1.0  |1.0
|23 |1.0  |1.0  |1.0
|25 |1.0  |1.0  |1.0
|27 |1.0  |1.0  |1.0

La siguiente tabla muestra las tasas de acierto para 10 segmentos con duracion de segmentos t = 4.5, 7, 12, 15, 17, 20, 23, 25 y 27 correspondientes.


|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5|1.0  |1.0  |1.0
|7  |1.0  |1.0  |1.0
|12 |1.0  |1.0  |1.0
|15 |1.0  |1.0  |1.0
|17 |0.9  |1.0  |1.0
|20 |1.0  |1.0  |1.0
|23 |1.0  |1.0  |1.0
|25 |1.0  |1.0  |1.0
|27 |1.0  |1.0  |1.0


La siguiente tabla muestra las tasas de acierto para 20 segmentos con duracion de segmentos t = 4.5, 7, 12, 15, 17, 20, 23, 25 y 27 correspondientes.


|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5|0.9  |1.0  |1.0
|7  |0.95 |1.0  |1.0
|12 |1.0  |1.0  |1.0
|15 |1.0  |1.0  |1.0
|17 |0.9  |1.0  |1.0
|20 |1.0  |1.0  |1.0
|23 |0.95 |1.0  |1.0
|25 |1.0  |1.0  |1.0
|27 |1.0  |1.0  |1.0


con 30 segmentos:

|T[seg]| 0dB |10dB|20dB|
|:-:   |:--: |:--:|:--:|
|4,5 |0.833333  |1.0  |0.966667
|7   |0.833333  |1.0  |1.0
|12  |0.966667  |1.0  |1.0
|15  |1.0       |1.0  |1.0
|17  |1.0       |1.0  |1.0
|20  |1.0       |1.0  |1.0
|23  |1.0       |1.0  |1.0
|25  |1.0       |1.0  |1.0
|27  |0.966667  |1.0  |1.0
=#


# ╔═╡ 6d76f2f2-3e67-11eb-04dc-0580a2072dda
md"""
#### Ejercicio 14) *(OPTATIVO)*

**Reproduzca y grabe un segmento de alguna de las canciones e intente reconocerla con su algoritmo. Puede utilizar dispositivos como el micrófono de su PC, su teléfono celular, una radio, etc. Comente los resultados obtenidos así como las condiciones de ruido de fondo y los dispositivos utilizados. (Recuerde que su función debe recibir audios a 44100 Hz.)**
"""

# ╔═╡ e8200592-3e7a-11eb-0711-ddf863314bee


# ╔═╡ 685698fa-3e67-11eb-2698-937dd4801b5c
md"""
#### Ejercicio 15) *(OPTATIVO, SOLO PARA LOS MÁS ENTUSIASTAS)*

**Repita el ejercicio 12 para T=10 solamente, afectando los segmentos con otro tipos de distorsiones elegidas entre las siguientes:**

* **Saturación**
* **Cuantización**
* **Cambio de velocidad**
* **Ecualización**

**Si elige saturación y/o cuantización, muestre el efecto que ocasiona sobre el audio mediante un espectrograma. Justifique si estas distorsiones pueden considerarse o no sistemas LTI.**
"""

# ╔═╡ ea11dfe2-3e7a-11eb-19fa-db2ebfcfecdc


# ╔═╡ 62b03a84-3e67-11eb-3949-2dc573c7d956
md"""
#### Ejercicio 16) *(OPTATIVO; SÓLO PARA LOS MÁS OBSESIVOS)*

**Verifique cómo cambia la tasa de aciertos cuando:**

* **cambia el solapamiento**
* **incrementa el solapamiento sólo al momento de identificar pero no al armar la base de datos**
* **cambia la longitud de la ventana manteniendo la tasa de frames por segundo**
* **cambia el algoritmo de extracción de características por algún otro que se le ocurra**
"""

# ╔═╡ 09446236-3e7b-11eb-2872-d720f72de7ee


# ╔═╡ 562997ce-3e67-11eb-015a-d318429ed230
md"""
# Apéndice

## La estructura de la base de datos

La base de datos es un vector de $2^{20}$ elementos, cada uno de los cuales es un vector de tamaño dinámico. Para cada frame de la huella acústica, se genera un valor a guardar en estos vectores de tamaño dinámico. El índice en el que se guarda cada valor se obtiene pasando a decimal el número binario conformado por las características del frame en cuestión. Note que la cantidad de elementos del vector de la base, $2^{20}$ se debe a que cada frame de la huella acústica posee 20 elementos binarios.

Cada valor a guardar es un objeto de tipo `FrameID`, que consiste de un entero de 32 bits sin signo, que determina el número de frame, y un entero de 16 bits sin signo que guarda el índice de la canción. Observe que el máximo número de canciones distintas que se podrán almacenar en esta tabla es $2^{16}$. Obserbe también que múltiples FrameIDs de distintos frames pueden requerir ser guardados en la misma posición por tener la misma huella; por esto, cada elemento es un vector que puede crecer según se requiera.

Esta implementación puede ser optimizada pero es un balance razonable entre eficiencia y simpleza, suficiente para este trabajo.


## La búsqueda de la canción más probable

La búsqueda en la tabla *hash* se realiza partiendo de una huella acústica de entrada y devuelve el ID de la canción más probable a la cual corresponde como salida. 

Para cada *frame* de la huella se obtiene el índice del vector de la base de datos que debe consultarse, pasando la característica de cada *frame* de binario a decimal, del modo inverso que cuando se almacenan elementos en la base, y se extraen todos los elementos almacenados en esa posición. Cada uno de estos elementos posee el ID de la canción a la cual corresponde y el número del frame que le correspondía originalmente.

Una forma sencilla de decidir a qué canción corresponde la huella podría ser elegir, de todos los elementos que se extrajeron, el ID que aparezca el mayor número de veces.

 Un refinamiento al criterio anterior es, dado que para cada ID se dispone del número de *frame* del cual se extrajo, quedarse con la mayor cantidad de ID dentro de un intervalo de *frames* igual a la longitud de la huella que se está consultando. Esto es lo que está implementado en la función *query_DB*. 

Por ejemplo, supongamos que se realiza un query a la base de datos con una huella que posee una longitud de 5 frames y se obtienen los siguientes 7 elementos:


![sys-tpe7.png](https://cdn.nextjournal.com/data/QmWrqDXaxxMXmbQ56JmTwH8EFfAPYmkB6uGQRwtqgBCXtH?filename=sys-tpe7.png&content-type=image/png)

Bajo el primer criterio se declararía ganador al ID 7 dado que aparece mayor cantidad de veces, pero bajo el segundo criterio se decide por el ID 3, debido a que en el intervalo de 5 frames que tiene la huella de entrada el ID 7 aparece como máximo 2 veces y el ID 3 aparece 3 veces.

# Bibliografía

![sys-tpe8.png](https://cdn.nextjournal.com/data/QmcYW98tYgXiX9gm4dzpWvQEjBW8dbr5bnwHkdUwRWvZE2?filename=sys-tpe8.png&content-type=image/png)
"""

# ╔═╡ Cell order:
# ╟─09062294-3e5f-11eb-176f-dfcbf841f111
# ╟─8f1394ee-3e63-11eb-0093-e75468460dc5
# ╟─7c04611e-3e61-11eb-0aa5-eb97132ace53
# ╟─82bdb9f0-7a26-11eb-0df3-0553dd890af1
# ╟─adc46380-3e63-11eb-2422-5bfe1b5052ba
# ╠═a3bf22c4-3ea3-11eb-3d3d-adfdfc171c33
# ╟─54e21790-7a14-11eb-14f3-9bf2cbf17fe6
# ╠═d132a762-3ea3-11eb-3494-692576a31f34
# ╠═aeea5400-794c-11eb-35b4-b763196181af
# ╟─4612d502-7961-11eb-0a80-11d812178ad8
# ╠═28c5ed26-3e6b-11eb-1d44-01e209b92f00
# ╟─92649f90-73cd-11eb-0df8-4958d753607d
# ╠═a2fa88b0-73cd-11eb-1336-9fbf72b0ddd8
# ╟─b65f8000-7962-11eb-1ba6-213e0aac1846
# ╟─a7727c3e-73cd-11eb-3bdc-8dd64b6b43ad
# ╟─a83a1200-73cd-11eb-353c-751703d316cd
# ╟─a56414d0-7621-11eb-0140-e57ff8c679ff
# ╟─1e7678f0-73ce-11eb-02a6-053cc2e84b2f
# ╠═202f4c80-73ce-11eb-2b03-1dc6c7240993
# ╟─b9ad22ac-3e67-11eb-35e1-7f4579b64838
# ╟─2b430db0-7963-11eb-3585-893f8168af1c
# ╠═4fc8c804-3ea4-11eb-3e97-eb6709f1c0aa
# ╟─56f4c0b0-73ce-11eb-09fd-fba3b409e4d9
# ╟─d987dc10-7621-11eb-26d9-d148f796fc11
# ╟─7776221e-73ce-11eb-003c-857b22a7e7c1
# ╠═bffcb902-73ce-11eb-31c9-05d7c512a3d0
# ╟─b60ae59e-3e67-11eb-123e-11c0cba7d09e
# ╟─de8b8c70-7964-11eb-3127-d3731697c053
# ╠═4e904a84-3ea4-11eb-0c12-b1fccd5f7036
# ╟─98aab180-73cf-11eb-092d-cb520eba3c36
# ╠═d61f59b0-73d1-11eb-0edb-938b1f3defa2
# ╟─f538d140-73d2-11eb-2248-abda76c139a5
# ╟─c0621930-796e-11eb-1a94-a118feca6e6b
# ╟─1ae45790-7a1b-11eb-03b9-a3212ee88756
# ╟─b2025250-3e67-11eb-39a2-73292bbf17c9
# ╟─4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
# ╟─3936fb00-73ed-11eb-25cd-812cdcf591fd
# ╠═c3a33cae-73d7-11eb-2707-d3604f680c0c
# ╟─1500ddb0-73e7-11eb-024b-1f68d5ac7525
# ╟─b07a0860-79d3-11eb-197f-4ba2019e2521
# ╟─5fef73e0-73e7-11eb-36e6-d1a47f2c456c
# ╠═27b2c6d0-73d8-11eb-33d5-e779bb4bf27d
# ╠═79411da0-79ee-11eb-3f3d-5b00dcffd43c
# ╠═ece5f452-79ef-11eb-0b5f-8141791cbab8
# ╠═54934580-73d8-11eb-2211-3f5357c88358
# ╠═340d0210-73d8-11eb-323c-e9501bdf682f
# ╟─fbde8120-79ef-11eb-05f7-2def33cc8263
# ╠═7e90a16e-73d8-11eb-1e23-0b1e366467aa
# ╟─c4921f20-79ee-11eb-1cc5-99c64e998496
# ╠═b18648f0-73d8-11eb-2f64-ff6fb23b096a
# ╠═ff2934c0-79ee-11eb-05ea-b77dfd7f1a73
# ╠═92bed250-75e2-11eb-3e64-d15f6a219026
# ╟─9fee5cf0-75f8-11eb-0404-cd8ad839260c
# ╠═2fcff8a0-761d-11eb-12d0-1736a6b45e85
# ╟─29aa4a30-79f0-11eb-207d-6f7e985947d0
# ╠═6ef5b922-73e3-11eb-0504-652233cd6dda
# ╟─d08961ee-7550-11eb-0529-513ae044a205
# ╠═72acf4c0-73e3-11eb-1975-77a7c0a9c562
# ╟─c3f1fca0-79f0-11eb-287e-4909d339ae38
# ╠═7c44b4a0-73e3-11eb-03c4-afb2354639a0
# ╟─cc11a150-73e3-11eb-0459-4fedbabe462d
# ╠═16be5cc0-7777-11eb-1b8e-b765ed6d37cd
# ╟─af4f3da4-3e67-11eb-3cc6-3378e0c12667
# ╠═3aa5434e-3ea4-11eb-20aa-b15564d4eb90
# ╟─38c05a80-73e4-11eb-274d-d538d2e3fb65
# ╟─982538c4-3e67-11eb-229e-dd2531a540d6
# ╠═f3971b40-73ea-11eb-3f2e-dd3a4443b9ab
# ╠═f906add0-74b1-11eb-156d-9b7c54bfd2ba
# ╠═39f5fc86-3ea4-11eb-37f3-25feb7d2aee6
# ╟─f7de1ba0-79fc-11eb-1fa8-c9f36f88d7ba
# ╟─e1906f20-7a00-11eb-3204-2f8cac355516
# ╠═5cc77020-74bc-11eb-154a-277125d7a831
# ╟─fd9d3860-7a00-11eb-2d95-b5700362c015
# ╟─aa6bbab2-7620-11eb-36d0-9b82996620bb
# ╟─987a7410-7a01-11eb-1f48-dd5f50a1cd81
# ╟─b52628a0-7620-11eb-1dda-65bc93222c5c
# ╠═6d1698d0-74a7-11eb-15bf-35d2eae0061f
# ╟─9309e284-3e67-11eb-1ab2-612f6c748c3b
# ╟─e7793b50-7a01-11eb-35ce-bfc0ff756a76
# ╠═5f636b02-3ea4-11eb-3f78-6f693a936992
# ╟─59fddd50-7513-11eb-1f17-07ef0c9f45f0
# ╠═6813ab40-7513-11eb-11b4-8f3040b56c91
# ╟─22d04590-7a02-11eb-3fe2-ab6b63005994
# ╠═a46234e0-7513-11eb-21d2-8592c8fec31d
# ╠═aeb41d00-7513-11eb-0980-9b6250c020af
# ╟─61106600-7a02-11eb-256f-d724de274557
# ╠═e6a127a0-751b-11eb-34ac-418d1ad01639
# ╟─d3fe3e30-7a02-11eb-25ea-85d66902a047
# ╠═46b25e4e-7609-11eb-3ef8-bff9a5a5898d
# ╟─7fda0ee0-752c-11eb-2b7d-5bba9b2a3592
# ╟─8a696ed0-7524-11eb-0764-214ebea3d1e7
# ╟─0fa90dc0-7a03-11eb-3301-5166ff8dfc61
# ╠═4a342420-75df-11eb-35be-0bc9cc550163
# ╟─8deaf928-3e67-11eb-0327-31e0f74de814
# ╟─198f8280-7a06-11eb-2a86-d16bc498aaef
# ╠═6476e9fc-3ea4-11eb-3873-b765108f4bab
# ╟─2435f260-7618-11eb-1b22-d3eda517d5f0
# ╠═d33d59f0-75f6-11eb-2c74-1b8e3be2f9c1
# ╟─feb5ec9e-7a06-11eb-034a-c169e0377755
# ╠═5bf232a0-7568-11eb-302c-a17935859a8a
# ╠═42c730b0-7625-11eb-1c4e-6b516c1eb2e7
# ╟─7b0ea820-7771-11eb-2c41-15bf560369ab
# ╟─89743a62-3e67-11eb-209e-9b1f3cc84e34
# ╟─d4bbda22-7797-11eb-1835-ddf37a08cc0f
# ╟─304f1ef0-7776-11eb-1926-13a4a3e50e26
# ╠═76604f60-777e-11eb-2543-09be49dcc4b1
# ╟─bd363c20-7a07-11eb-1a07-a1e3ea91852f
# ╠═d4179090-77a7-11eb-3c1b-4da6e8152c16
# ╠═4b0b7360-7776-11eb-1bd5-ddd0c82e9e06
# ╟─5795e360-7a08-11eb-0142-31c7502d8428
# ╠═741304fe-3ea4-11eb-15e8-09908d98ecb3
# ╠═ed9dce70-7785-11eb-27b2-8d22fa913c96
# ╟─855a7d2e-3e67-11eb-0f46-a5c786d5caf3
# ╟─4df7d760-7796-11eb-2563-932eb8c9fece
# ╠═7a02bdf0-7793-11eb-3606-9168a11719a4
# ╠═c19daa60-7795-11eb-03a6-d134602f316a
# ╠═e66049e2-7793-11eb-0560-1189121d97f9
# ╟─81717fc8-3e67-11eb-05fc-5bde46597f8a
# ╟─8a2d6fe0-7a09-11eb-1946-694ee29e9d86
# ╠═b91537ac-3ea4-11eb-14d6-d341c535d83e
# ╟─73333e92-3e85-11eb-26b6-7f0309ef2ee9
# ╟─feb5d512-3e85-11eb-0116-29e4d9539595
# ╠═0cf7ba9c-3e74-11eb-18e2-c38aa20f9e9a
# ╠═928af020-7a0b-11eb-325c-81c84be31414
# ╠═e9255b8c-3e74-11eb-2960-5d01b0c99b13
# ╠═5d0a5800-779a-11eb-2f2c-335d6f94ef06
# ╟─7c7c1424-3e67-11eb-1da0-5dbad0171b20
# ╠═415e32e6-3e76-11eb-17fa-23bd653fb975
# ╟─76ce23dc-3e67-11eb-0be0-91b6781840fb
# ╟─a9b2ac20-7a0b-11eb-13cd-2de6007fc26d
# ╠═7d82ea60-779b-11eb-2e75-8997f2c88b53
# ╟─c970df70-77ac-11eb-384a-9b8a5e2e9547
# ╠═b99bea60-7601-11eb-3f06-ef522296b39d
# ╠═00da49d0-7602-11eb-31f1-dfe5b07d6349
# ╟─7229577a-3e67-11eb-0c71-f383056175d1
# ╟─746ccbb0-7a0e-11eb-11f0-f107a2b47442
# ╠═841f7430-77b2-11eb-180b-11a3ec034424
# ╠═8fe96ca0-77ce-11eb-1d58-27201d27fc2f
# ╠═129267b0-77b6-11eb-2fa3-5914fd885d43
# ╟─7b6294de-77bb-11eb-2111-01df80d418e5
# ╟─713783c0-77bd-11eb-2b34-17a3651f41c9
# ╠═612586c0-7a25-11eb-247f-8f55aa84cb91
# ╟─8274ebc0-7a13-11eb-2632-3f936f3eb7e2
# ╠═6d76f2f2-3e67-11eb-04dc-0580a2072dda
# ╠═e8200592-3e7a-11eb-0711-ddf863314bee
# ╟─685698fa-3e67-11eb-2698-937dd4801b5c
# ╠═ea11dfe2-3e7a-11eb-19fa-db2ebfcfecdc
# ╟─62b03a84-3e67-11eb-3949-2dc573c7d956
# ╠═09446236-3e7b-11eb-2872-d720f72de7ee
# ╟─562997ce-3e67-11eb-015a-d318429ed230
