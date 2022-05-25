import numpy as np

def bin2gray(bi):
    
    order = bi.shape[1]
    lin = bi.shape[0]
    g = np.zeros((lin, order))

    idx = 0
    for bits in bi:
        g[idx, 0] = int(bits[0])
        
        for l in range(1, order):
            g[idx, l] = int(bits[l] ^ bits[l - 1])
        
        idx = idx + 1

    return g

def qpskmodulator(input_bits, average_power=1, phase_offset=np.pi/4):

    symbols = np.array([1, 1j, -1, -1j]) * np.exp(1j * phase_offset)
    symbols = np.sqrt(average_power) * symbols

    input_arr = input_bits.reshape((int(len(input_bits) / 2)), 2)
    gray_code = bin2gray(input_arr)

    output = []
    for symbol in gray_code:
        s = int(2 * symbol[0] + symbol[1])
        output.append(symbols[s])

    return np.array(output)

def awgn(sig, SNR, sig_power=1):

    SNR_linear = 10**(SNR / 10)
    noise_power = sig_power / SNR_linear

    dim = sig.shape
    if not np.isreal(sig[1, 1]):
        noise = np.sqrt(noise_power / 2) * (np.random.normal(0, 1, dim) +
                                            1j * np.random.normal(0, 1, dim))
    else:
        noise = np.sqrt(noise_power) * np.random.normal(0, 1, dim)

    sig_noisy = sig + noise

    return sig_noisy