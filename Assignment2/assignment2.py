import numpy as np;
from scipy.signal import zpk2tf,freqz;
import matplotlib.pyplot as plt;
from scipy import signal;
from scipy.io import wavfile as wav;
from scipy.fft import fft,fftfreq;
from numpy import log10;
import warnings
warnings.filterwarnings("ignore")

def pre_emphasize(sound,alpha):

    pre_sound = np.zeros(len(sound))
    for i in range(1,len(sound)):
        pre_sound[i-1] = sound[i] - alpha*sound[i-1]
    return pre_sound


def plot_sound(sound,title_):
    plt.figure()
    plt.plot(sound)
    plt.title(title_)
    plt.grid()
    plt.xlabel('Samples')
    plt.ylabel('Magnitude')
    plt.savefig('outputs/'+title_+".png")


def plot_spectrum(dft_,freq,names,title_,plots):

    colors = ['b','g','k','r']
    plt.figure()
    for i in range(plots):
        len_ = int(len(dft_[i])/2)                
        plt.plot(abs(freq[i]),20*log10(abs(dft_[i])),colors[i],linewidth=0.7)
    
    plt.legend([i for i in names])
    plt.xlim(xmin=0)
    plt.grid("True")
    plt.ylabel("Magnitude(dB)")
    plt.xlabel("Samples(n)")
    plt.title(title_)
    plt.savefig('outputs/'+str(plots)+title_+".png")


def LPRecursion(window_signal,p):
    r= np.correlate(window_signal_pre,window_signal_pre,mode='full')
    r= r[-len(window_signal_pre):]  # keeping only the positive coefficients since correlation is symmetric

    e = np.zeros(p+1)
    a = np.zeros((p+1,p+1))
    g = np.zeros(p+1)

    e[0] = r[0]
    g[0] = np.sqrt(e[0])
    a[1][0] = 1

    k=r[1]/e[0]
    a[1][1] = k
    e[1]=(1-k**2)*e[0]
    g[1] = np.sqrt(e[1])

    for i in range(2,p+1):    
        temp = 0
        for j in range(1,i):
            temp+=a[i-1][j] *r[i-j]
        k = (r[i] - temp)/e[i-1]    
        a[i][i] = k    
        for j in range(1,i):
            a[i][j] = a[i-1][j] - k*a[i-1][i-j]    
        e[i] = (1-k**2)*e[i-1]
        g[i] = np.sqrt(e[i])
        a[i][0]=1
    
    return g,e,a

def poleZeroPlot(b,a,title_):

    poles = np.zeros_like(a)
    poles[0]=1
    poles[1:len(a)] = -a[1:len(a)];
    b=[b]
    z, p, k = signal.tf2zpk(b, poles)
    p = p[p!=0]
    
    fig = plt.figure(figsize=(5,5))
    ax=fig.add_subplot(1, 1, 1)
    plt.title(title_)
    plt.plot(np.real(z), np.imag(z), 'ob')
    plt.plot(np.real(p), np.imag(p), 'sr',markersize=5,fillstyle="full")
    circ = plt.Circle((0, 0), radius=1,facecolor='None',color='black', ls='solid', alpha=0.1)
    ax.add_patch(circ)
    plt.axhline(0,color='black',alpha=0.4)
    plt.axvline(0,color='black',alpha=0.4)
    plt.ylim((-2.0, 2.0))
    plt.xlim((-2.0,2.0))
    plt.legend(['Zeros', 'Poles'])
    plt.ylabel('Real')
    plt.xlabel('Imaginary')
    plt.grid()

    plt.savefig('outputs/'+str(title_))


def plot_LPC(p,a,g,dft_pre,fs,freq):
    for i in range(1,p+1):
        poles = np.zeros_like(a[i])
        poles[0] = a[i][0]
        poles[1:] = -a[i][1:len(a[i])+1] 
        w,h = freqz(g[i],poles)

        plt.figure()
        plt.title("LPC spectrum for poles="+str(i))
        plt.plot(abs(freq),20*log10(abs(dft_pre)),'g',linewidth=0.7)
        plt.plot(w*fs/(2*np.pi),20*log10(abs(h)),'r',linewidth=2)
        plt.ylabel("Magnitude(dB)")
        plt.xlabel("Frequency(w)")
        plt.xlim(xmin=0)
        plt.grid("True")
        plt.legend(['Hamming window output' , 'LP analysis estimate'])
        plt.savefig('outputs/LPCSpectrum_poles_'+str(i)+".png")

def plot_acf(acf,title_,fig_name):

    plt.figure()
    plt.plot(acf)
    plt.xlabel('Lag(Samples n)')
    plt.ylabel('Autocorrelation')
    plt.grid()
    plt.title(title_)
    plt.savefig('outputs/'+fig_name+'.png')



if __name__=='__main__':

    sound_file = wav.read('aa.wav')
    fs,sound = sound_file[0],sound_file[1]
    alpha=0.98
    
    pre_sound = pre_emphasize(sound,alpha)
    plot_sound(sound,'Original Sound Plot')
    plot_sound(pre_sound,'Pre-Emphasized Sound Plot')

    
    win_length = 30
    window_size = int(win_length*fs/1000)
    center = int(len(pre_sound)/2)
    window_signal_pre = pre_sound[center - window_size : center] * np.hamming(window_size)
    window_signal = sound[:window_size] * np.hamming(window_size)
    
    dft_orig = [fft(window_signal, n=1024)]
    dft_pre = [fft(window_signal_pre, n=1024)]
    freq = [fftfreq(dft_orig[0].shape[-1], 1/fs)]

    combined = [*dft_orig,*dft_pre]
    freq = [*freq,*freq]
    plot_spectrum(dft_orig,freq,["Original"],title_="Hamming Window magnitude response(Original Sound)",plots=1)
    plot_spectrum(dft_pre,freq,["Pre Emphasized"],title_="Hamming Window magnitude response(Pre-Emphasized Sound)",plots=1)
    plot_spectrum(combined,freq,["Original","Pre Emphasized"],title_="Combined output - Pre-emphasized and original",plots=2)

    
    p=10
    g,e,a = LPRecursion(window_signal_pre,p)

    plt.figure()
    plt.title("Error signal energy plot")
    plt.plot(20*log10(abs(e)),'b-*',linewidth=0.7)
    plt.ylabel("Energy error signal(dB)")
    plt.xlabel("Number of poles(p)")
    plt.xlim(xmin=0)
    plt.grid("True")
    plt.savefig('outputs/Error_signal.png')


    poleZeroPlot(g[6],a[6][:],'Pole-Zero plot for p=6')
    poleZeroPlot(g[10],a[10][:],'Pole-Zero plot for p=10')


    order =10

    plot_LPC(order,a,g,dft_pre[0],fs,freq[0])

    
    ###Reconstruction of the signal

    input_signal = window_signal_pre
    zeros_filter = a[order][:]
    poles = g[order]
    output = np.zeros(len(input_signal))

    for i in range(len(input_signal)):
        temp = 0
        for j in range(0,order):
            if i-j>=0:
                output[i]+=(-zeros_filter[j]*input_signal[i-j])
        output[i] = output[i]/poles

    acf = np.correlate(output,output,mode="full")
    original_acf =np.correlate(window_signal_pre,window_signal_pre,mode="full")
    plot_acf(acf,'Autocorrelation function for the estimated signal','acf_estimated')
    plot_acf(original_acf,'Autocorrelation function for the estimated signal','acf_original')



    #### Reconstructing the original signal form the estimated filter by inverse filtering
    estimated_f = 136
    t = np.linspace(0, 300, int(300*fs/1000), endpoint=False)
    sig = (signal.square(2 * np.pi *(estimated_f/1000)* t[0:1000],duty=0.08)+1)/2
    plt.figure()
    plt.plot(t[0:1000], sig)
    plt.title('Impulse train')
    plt.xlabel('Samples(n)')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.savefig('outputs/impulse_train.png')


    order =10
    output=np.zeros(sig.shape)

    for i in range(len(sig)):
        output[i] = g[order]*sig[i]
        for j in range(0,order):
            if i-j>=0:
                output[i]+=a[order][j]*output[i-j]

    output_pre = np.zeros(sig.shape)
    alpha=0.98
    for i in range(1,len(output)):
        output_pre[i] = output[i] + alpha*output_pre[i-1]
        
    maxi = np.max(output_pre)
    output_pre = output_pre/maxi

    plt.figure()
    plt.plot(output_pre)
    plt.title('Estimated output sound wave')
    plt.xlabel('Samples(n)')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.savefig('outputs/estimated_output.png')
    wav.write('outputs/result.wav',fs,output_pre)



    