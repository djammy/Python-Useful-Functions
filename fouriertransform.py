from scipy import fftpack

def FourierTransform(time,amplitude):
	Spectrum = fftpack.fft(amplitude)
	f_s = 1/(time[1]-time[0])
	Frequency = fftpack.fftfreq(len(time)) * f_s
	return(Frequency,Spectrum)