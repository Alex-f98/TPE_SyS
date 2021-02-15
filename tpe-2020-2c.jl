### A Pluto.jl notebook ###
# v0.12.20

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

# ╔═╡ 8f1394ee-3e63-11eb-0093-e75468460dc5
md"# Trabajo Práctico Especial de Señales y Sistemas"

# ╔═╡ 7c04611e-3e61-11eb-0aa5-eb97132ace53
html"""
<h2 class="pm-node nj-subtitle">Identificación de canciones mediante huellas digitales acústicas</h2>
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

# ╔═╡ d132a762-3ea3-11eb-3494-692576a31f34
function load_audio(fn)
	
    xraw = load(fn)
	
	# Chequeo que la frecuencia de muestreo sea "sr = 44100hz"
    @assert samplerate(xraw) == sr  samplerate(xraw) 

    return collect(xraw)
end

# ╔═╡ 3da199b0-6b3d-11eb-04e6-fff2ac06d553
#toma un audio (.ogg) de dos canales, los promedia a un solo canal.
function to_mono(fn) 
  
  stereo = load_audio(fn)
  mono = vec( mean(stereo, dims=2) )
  
  return mono #collect(mono)
end

# ╔═╡ 53a1ecb0-6b3d-11eb-1f24-55b65316c313
x = to_mono("Pink.ogg")   #..maldita vida...

# ╔═╡ 09062294-3e5f-11eb-176f-dfcbf841f111
begin
	import Pkg
	Pkg.add.(["DSP", "FFTW", "FileIO", "WAV", "PlutoUI", "Plotly"])
	Pkg.add.(["LibSndFile", "SampledSignals", "Statistics"])
	Pkg.add(["ImageMagick", "Images", "ImageIO", "Plots", "IterTools"])
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

# ╔═╡ 0a039420-6d75-11eb-11e6-3579dc60d9bb
md"""Notar que hay 441.000 muestras del tipo Float32 un con "sr=44100" mi frecuencia de muestreo el video debe durar 10 segundos"""

# ╔═╡ 28c5ed26-3e6b-11eb-1d44-01e209b92f00
# Wrapper para ver un vector de muestras de audio a 44100 Hz como un widget escuchable
sound(x) = SampleBuf(x, sr);

# ╔═╡ 621134d0-6d78-11eb-2959-99d0c3971552
sound(x)

# ╔═╡ 8e98a8ee-6d76-11eb-24d0-5305b72b633c
md""" Se grafica parte de la señal de audio se toma 1500 muestras de las 441.0000 """ 

# ╔═╡ 117bab80-6d79-11eb-235b-67753330b843
n = 1:1:1500;

# ╔═╡ 7f346530-6d7a-11eb-1066-75652a808a15
tseg = n./44100;

# ╔═╡ 344bedb0-6d73-11eb-3ec5-158ce2740b64
plot(tseg,x[n],
	label =  "f = 44100",
	legend = true,
	xlabel = "Tiempo(Seg)",
	ylabel = "Amplitud",
	#xlims=(0, 1500),
	#ylims=(-1.1, 1.1)

)

# ╔═╡ b9ad22ac-3e67-11eb-35e1-7f4579b64838
md"""
#### Ejercicio 2)

**Mediante la función `fft`, realice un análisis espectral de la señal mostrando la densidad espectral de potencia (PSD) en función de la frecuencia en Hz.  Verifique que la mayor PSD se concentra en las bajas frecuencias. Utilice escala logarítmica para mostrar esultado.**

**La PSD, ``S_{XX}``, la puede estimar como ``S_{X X}(\omega) = 1/T |X(\omega)|^2``, donde ``T`` es la duración del segmento temporal que usa para la estimación, o mejor aún, dividiendo al audio en segmentos, realizando la estimación anterior con cada uno, y finalmente promediándolas.**
"""

# ╔═╡ 4fc8c804-3ea4-11eb-3e97-eb6709f1c0aa
function PDS(x)
	ft = fft(x)
	sxx = 1/(t.*(	abs.(ft)	))
end

# ╔═╡ 3f5d6960-6d80-11eb-1eab-c53c17a520ef

paso=length(x)/100;

# ╔═╡ 9450ce92-6d7f-11eb-259c-ddf6325b6d7f
tst = x#c = [ fft(x[i:i+paso]) for i in 1:paso:length(x) ];

# ╔═╡ ee59e260-6d88-11eb-3256-c36cab60d6b0
#tener en cuenta que fft te da la transformada discreta de fourier, yo quiero la continua. por la relacion que hay yo tengo que dividir por la frecuencia de muestreo.
## retorna el PSD en dB.
psd = let	
	n =	1000;		#4410 						#cantidad de muestras por segmento
	t_win = n / sr;     					    #ancho del segmento.
	#fft(tst_seg) ./ sr == X(W).
	mean(10log10.(abs.( fft(tst_seg) ./ sr ).^2 ./t_win)
	for tst_seg in arraysplit(tst, n, 0)
	)
#arraysplit me devuelve algo asi como un "vector de vectores de largo n".
end

# ╔═╡ b09f0110-6d94-11eb-084a-6129c8fca2fa
#f_psd1 =range(1, sr/2, length =div(length(psd),2))

# ╔═╡ d0993fe2-6d8e-11eb-02af-f32e42bfe853
f_psd = let
	m = 1:length(psd)/2;
	(m.*sr)./length(psd) 
end

# ╔═╡ 211f97de-6d97-11eb-0408-85c92eeb0a0e
md""" preguntar si esto esta bien!!!!!!! EL EJE DE FRECUENCIA(f_psd) NO ME DA MUCHA CONFIANZA"""

# ╔═╡ 84d3ce80-6d8a-11eb-09f9-81d1392a1051
plot(f_psd,
	psd[1:div(length(psd),2)],
	xlims = (0, 6000),
	ylims =	(-80, -40),
	title =  "PSD vs Frec.",
	xlabel = "Frecuencia (Hz)",
	ylabel = "PSD en dB (1/Hz)",
	legend = false
	)

# ╔═╡ b60ae59e-3e67-11eb-123e-11c0cba7d09e
md"""
#### Ejercicio 3)

**Mediante la función `fft`, obtenga el espectro de una ventana rectangular y una de Hamming de igual duración y grafique su potencia en escala logarítmica en un mismo gráfico para compararlos. ¿Cuál es la resolución en frecuencia de cada ventana? Al realizar un espectrogama, ¿qué ventaja tiene utilizar la ventana de Hamming frente a la rectangular?**
"""

# ╔═╡ 0ae0430e-6d76-11eb-34c4-7b4ab551f242
let
  normalize(x) = x ./ sum(x)
   plot()
   plotwin(win) = plot!(log.(abs.(fft(padright(normalize(win(41)), 501))));
	legend= true,
   	xlims=(0, 150), ylims=(-15, 2), ylabel="dB")
  
  foreach(plotwin, [ones, hamming],)
  plot!()

end

# ╔═╡ 4e904a84-3ea4-11eb-0c12-b1fccd5f7036
let
	n = 41;
	sh= 501;
	
	normalize(x) = x ./ sum(x);
	log_fft(win) =	 plot(log.(abs.(fft(padright(normalize(win(n)), sh)))));
	
	plot(log_fft(ones), log_fft(hamming),
		label= ["ventana" "hamming"],
		legend= true,
   		xlims=(0, 250), ylims=(-15, 2), 
		ylabel="Amplitud (dB)"
	)
#es necesiario pasar a omega?...omega = (2.pi.k)/N.
end

# ╔═╡ eca8b002-6f35-11eb-2ef7-c1d8fb43e4cb
specplot(tst)

# ╔═╡ 4db2d110-6f30-11eb-187f-d36624b5dbfa
md""" Se puede observar la diferencia entre los lobulos secundarios entre la ventana ideal y la de hamming, la relacion entre el primer lobulo y el segundo es mucho mejor con hamming."""

# ╔═╡ b2025250-3e67-11eb-39a2-73292bbf17c9
md"""
#### Ejercicio 4)
**Implemente un sistema para reducir la frecuencia de muestreo del audio, de 44100 Hz a 5512.5 Hz. Muestre un diagrama en bloques y justifique su diseño. Diseñe el filtro pasabajos mediante el método de ventaneo explicando y justificando las decisiones, de forma tal que su retardo sea menor a 1 ms. Graficar un diagrama de polos y ceros, respuesta en frecuencia (módulo y fase), respuesta al impulso, y retardo de grupo.**
"""

# ╔═╡ 4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
#plot(10log10.(fftshift(abs.(fft(tst)))))

# ╔═╡ 288dc2a0-6fa8-11eb-0f59-81c9f9871214
md""" La velocidad de muestreo puede ser reducida por un factor entero M= 8, haciendo una decimacion de la señal discreta generando una nueva secuencia $x_d[n]$.
			$x_d[n] = x_d[m * M]= x_c[n*MT]$,
con T:periodo de muestreo.

En el espectro esto se puede ver con amplitud 1/N y ancho de banda desde $-M.\Omega_m$   a $-M.\Omega_m$ entonces la reduccion en frecuencia de muestreo puede generar un traslape para $\Omega_N > \pi/M$.

Por ello para no perder informacion se realiza un filtro anti-aliasin antes del submuestreo.
"""

# ╔═╡ b09a2b70-6fcb-11eb-344e-abf9e54cb366
let
	img_path = "decimador.png"
	img = load(img_path)
end

# ╔═╡ 596ad140-6fd2-11eb-11ba-b542ccbd645e


# ╔═╡ d3f093f0-6fcd-11eb-0683-2352d0d55bbc
begin
	PI= 3.14159265;
	fc= 2756; 					#frecuencia de corte(hz) fc < fs/2 = 44100/2.
	OmegaC = (2*PI*fc)/sr;
end

# ╔═╡ 5b221e30-6fd2-11eb-054b-518034184d4a
begin
	H_ideal(Ω) = u(Ω + OmegaC) - u(Ω - OmegaC)
	plot(H_ideal, -π, π)

end

# ╔═╡ 6bcd16f0-6fc7-11eb-3e01-6747e7351a51
h_ideal(n, OmegaC) = (OmegaC/PI) * sinc(n *(OmegaC/PI))	#Omega = 2pi/fm

# ╔═╡ 95d5127e-6fcd-11eb-2770-3b888691f1ec
let ns = -40:40	
  	stem(ns, h_ideal.(ns,OmegaC))
	
end

# ╔═╡ e752c6d0-6fa9-11eb-2ebd-19c22fe34d5c
function submuestrear(x,M)
  y = zeros(length(x)*M) 					#me creo un vector para rellenar
  [y[1 + M*(i-1)]= x[i] for i=1:length(x)] 	# aca hago inline loops :)
  return y
end

# ╔═╡ acee07a0-6fb0-11eb-2a65-5d76cb22034b
#fvtol     t=muestras/FS

# ╔═╡ 3f345a30-6fc8-11eb-2a37-95cb5d649423
begin
	order = 231;
	ns    = range(-(order - 1) / 2; length=order);
	winds = hamming(order);
	hs    = h_ideal.(ns,OmegaC) .* winds;
	
	stem(0:order-1, hs)
end

# ╔═╡ 2ea03a40-6fd2-11eb-3a3a-1f34c6888318
let
  samples = 500
  hfs = fft(padright(hs, samples))
  fs = range(0; step=2pi/samples, length=samples)
  
  plot(fs, abs.(hfs); xlims=(0, π))
  plot!(H_ideal, 0, π)
end

# ╔═╡ 9e1b3d62-6fd3-11eb-3738-17d34b356dfe
#aca reduciomos tst en sr/8.
#tst_8 = let nss = 1:length(tst)
#	p1 = tst .* h_ideal.(nss,OmegaC) .* winds;
#	submuestrear(p1.* winds,8)
	
#end

# ╔═╡ af4f3da4-3e67-11eb-3cc6-3378e0c12667
md"""
#### Ejercicio 5)

**Realice un espectrograma de la señal original y la filtrada y verifique los efectos del filtrado. Indique el tipo y longitud de ventana utilizada.**
"""

# ╔═╡ 982538c4-3e67-11eb-229e-dd2531a540d6
md"""
### Espectrograma

Los parámetros con los que se realice el espectrograma tienen una influencia directa sobre la extracción de las características. La elección de la ventana (tipo y longitud) es el parámetro principal, ya que debería lograr una resolución razonable tanto en tiempo como en frecuencia que permita distinguir las características espectrales y temporales con claridad.

#### Ejercicio 6)

**Realice un espectrograma de la señal sub-muestreada y muestre una imagen (con zoom) de alguna región donde se vean las características espectrales de la señal dada la ventana utilizada. Muestre y justifique cómo cambia la visualización de las características con 3 diferentes longitudes de ventana y comente qué longitud utilizará en el algoritmo. Realice las comparaciones con ventana rectangular y de Hamming.**
"""

# ╔═╡ 39f5fc86-3ea4-11eb-37f3-25feb7d2aee6


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

# ╔═╡ 5f636b02-3ea4-11eb-3f78-6f693a936992


# ╔═╡ 8deaf928-3e67-11eb-0327-31e0f74de814
md"""
Finalmente, las características se obtienen mediante una función que opera sobre E(m, n) según:


![sys-tpe4.png](https://i.imgur.com/dGoWXWE.png)

En palabras, $H[m, n]$ es 1 si la diferencia de energía entre las bandas $m$ y $m+1$ para el *frame* actual n es mayor a la del *frame* anterior $n-1$. Experimentalmente, se verficó que estas características son robustas ante varios tipos de procesamientos y distorsiones del audio\[5\].

#### Ejercicio 8)

**Implemente el algoritmo de extracción de características, calculando la huella. Debido al efecto borde, $H[m, n]$ debería resultar una matriz de 20 filas. Ejecute el algoritmo sobre un segmento de audio y muestre una imagen de la huella digital acústica obtenida.**
"""

# ╔═╡ 6476e9fc-3ea4-11eb-3873-b765108f4bab


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

# ╔═╡ 741304fe-3ea4-11eb-15e8-09908d98ecb3


# ╔═╡ 855a7d2e-3e67-11eb-0f46-a5c786d5caf3
md"""

#### Ejercicio 10)

**Observe que la cantidad de elementos a guardar en la base de datos se incrementa conforme la longitud de las ventanas del espectrograma inicial disminuye, o el solapamiento entre ventanas se incrementa. Determine el solapamiento entre ventanas del espectrograma para obtener una densidad de aproximadamente 25 elementos por segundo y utilice este valor para el ejercicio siguiente.**
"""

# ╔═╡ 088ca198-3e74-11eb-0cf3-23a983165a0d


# ╔═╡ 81717fc8-3e67-11eb-05fc-5bde46597f8a
md"""

#### Ejercicio 11)

**Ejecute la función `generar_DB` para confeccionar la base de datos completa de su lista de canciones. Utilice al menos 40 canciones para llenar la base de datos. Puede usar la lista de canciones provista, y/o usar una lista de canciones propia. (Recuerde verificar que la frecuencia de muestreo de sus canciones sea de 44100 Hz).**
"""

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

# ╔═╡ e9255b8c-3e74-11eb-2960-5d01b0c99b13
# db = generar_DB(songs; dir=songsdir);

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

# ╔═╡ f91ed600-3ea4-11eb-3945-6999aaa4d0dd


# ╔═╡ 7229577a-3e67-11eb-0c71-f383056175d1
md"""
#### Ejercicio 13)

**Repita el ejercicio 12 sumando ruido a los segmentos de audio. Utilice la función `randn` para generar las muestras de ruido. Evalúe tasa de aciertos para $SNR =0dB$, $10dB$ y $20dB$, mostrando sus resultados en una tabla para 9 combinaciones de longitud temporal y ruido. Nota: $SNR=10 log_{10}(P_X/P_N)$ donde $P_X$ es la potencia media de la señal sin ruido, y $P_N$ es la potencia media del ruido sumado a la señal. Para el cálculo de la potencia media puede utilizar la función `var`, que estima la varianza de una señal, ya que las señales de audio no deberían componente continua o valor medio.**
"""

# ╔═╡ 4e913d8e-3ea6-11eb-25e6-e3d03de7b3e0


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
# ╟─adc46380-3e63-11eb-2422-5bfe1b5052ba
# ╠═a3bf22c4-3ea3-11eb-3d3d-adfdfc171c33
# ╠═d132a762-3ea3-11eb-3494-692576a31f34
# ╠═3da199b0-6b3d-11eb-04e6-fff2ac06d553
# ╠═53a1ecb0-6b3d-11eb-1f24-55b65316c313
# ╟─0a039420-6d75-11eb-11e6-3579dc60d9bb
# ╠═28c5ed26-3e6b-11eb-1d44-01e209b92f00
# ╠═621134d0-6d78-11eb-2959-99d0c3971552
# ╟─8e98a8ee-6d76-11eb-24d0-5305b72b633c
# ╠═117bab80-6d79-11eb-235b-67753330b843
# ╠═7f346530-6d7a-11eb-1066-75652a808a15
# ╠═344bedb0-6d73-11eb-3ec5-158ce2740b64
# ╟─b9ad22ac-3e67-11eb-35e1-7f4579b64838
# ╠═4fc8c804-3ea4-11eb-3e97-eb6709f1c0aa
# ╠═3f5d6960-6d80-11eb-1eab-c53c17a520ef
# ╠═9450ce92-6d7f-11eb-259c-ddf6325b6d7f
# ╠═ee59e260-6d88-11eb-3256-c36cab60d6b0
# ╠═b09f0110-6d94-11eb-084a-6129c8fca2fa
# ╠═d0993fe2-6d8e-11eb-02af-f32e42bfe853
# ╟─211f97de-6d97-11eb-0408-85c92eeb0a0e
# ╠═84d3ce80-6d8a-11eb-09f9-81d1392a1051
# ╟─b60ae59e-3e67-11eb-123e-11c0cba7d09e
# ╠═0ae0430e-6d76-11eb-34c4-7b4ab551f242
# ╠═4e904a84-3ea4-11eb-0c12-b1fccd5f7036
# ╠═eca8b002-6f35-11eb-2ef7-c1d8fb43e4cb
# ╟─4db2d110-6f30-11eb-187f-d36624b5dbfa
# ╟─b2025250-3e67-11eb-39a2-73292bbf17c9
# ╠═4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
# ╟─288dc2a0-6fa8-11eb-0f59-81c9f9871214
# ╟─b09a2b70-6fcb-11eb-344e-abf9e54cb366
# ╠═596ad140-6fd2-11eb-11ba-b542ccbd645e
# ╠═d3f093f0-6fcd-11eb-0683-2352d0d55bbc
# ╠═5b221e30-6fd2-11eb-054b-518034184d4a
# ╠═6bcd16f0-6fc7-11eb-3e01-6747e7351a51
# ╠═95d5127e-6fcd-11eb-2770-3b888691f1ec
# ╠═e752c6d0-6fa9-11eb-2ebd-19c22fe34d5c
# ╠═acee07a0-6fb0-11eb-2a65-5d76cb22034b
# ╠═3f345a30-6fc8-11eb-2a37-95cb5d649423
# ╠═2ea03a40-6fd2-11eb-3a3a-1f34c6888318
# ╠═9e1b3d62-6fd3-11eb-3738-17d34b356dfe
# ╟─af4f3da4-3e67-11eb-3cc6-3378e0c12667
# ╟─982538c4-3e67-11eb-229e-dd2531a540d6
# ╠═39f5fc86-3ea4-11eb-37f3-25feb7d2aee6
# ╟─9309e284-3e67-11eb-1ab2-612f6c748c3b
# ╠═5f636b02-3ea4-11eb-3f78-6f693a936992
# ╟─8deaf928-3e67-11eb-0327-31e0f74de814
# ╠═6476e9fc-3ea4-11eb-3873-b765108f4bab
# ╟─89743a62-3e67-11eb-209e-9b1f3cc84e34
# ╠═741304fe-3ea4-11eb-15e8-09908d98ecb3
# ╟─855a7d2e-3e67-11eb-0f46-a5c786d5caf3
# ╠═088ca198-3e74-11eb-0cf3-23a983165a0d
# ╟─81717fc8-3e67-11eb-05fc-5bde46597f8a
# ╠═b91537ac-3ea4-11eb-14d6-d341c535d83e
# ╟─73333e92-3e85-11eb-26b6-7f0309ef2ee9
# ╟─feb5d512-3e85-11eb-0116-29e4d9539595
# ╠═0cf7ba9c-3e74-11eb-18e2-c38aa20f9e9a
# ╠═e9255b8c-3e74-11eb-2960-5d01b0c99b13
# ╟─7c7c1424-3e67-11eb-1da0-5dbad0171b20
# ╠═415e32e6-3e76-11eb-17fa-23bd653fb975
# ╟─76ce23dc-3e67-11eb-0be0-91b6781840fb
# ╠═f91ed600-3ea4-11eb-3945-6999aaa4d0dd
# ╟─7229577a-3e67-11eb-0c71-f383056175d1
# ╠═4e913d8e-3ea6-11eb-25e6-e3d03de7b3e0
# ╟─6d76f2f2-3e67-11eb-04dc-0580a2072dda
# ╠═e8200592-3e7a-11eb-0711-ddf863314bee
# ╟─685698fa-3e67-11eb-2698-937dd4801b5c
# ╠═ea11dfe2-3e7a-11eb-19fa-db2ebfcfecdc
# ╟─62b03a84-3e67-11eb-3949-2dc573c7d956
# ╠═09446236-3e7b-11eb-2872-d720f72de7ee
# ╟─562997ce-3e67-11eb-015a-d318429ed230
