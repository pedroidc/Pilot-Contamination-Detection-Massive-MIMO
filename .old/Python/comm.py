import numpy as np


def qpskmodulator(input_bits, average_power=1, phase_offset=np.pi / 4):

    # Generate symbols
    symbols = np.array([1, 1j, -1, -1j]) * np.exp(1j * phase_offset)
    symbols = np.sqrt(average_power) * symbols

    input_arr = input_bits.reshape((int(len(input_bits) / 2)), 2)
    gray_code = bin2gray(input_arr)

    output = [
    ]  # Pre allocate in the future to prevent unnecessary memory usage
    for symbol in gray_code:
        #         print(symbol)
        s = int(2 * symbol[0] + symbol[1])
        output.append(symbols[s])

    return np.array(output)


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


def awgn(sig, SNR, sig_power=1):
    # SNR in dB

    SNR_linear = 10**(SNR / 10)  #SNR to linear scale
    noise_power = sig_power / SNR_linear
    # Noise power

    dim = sig.shape
    if not np.isreal(sig[1, 1]):
        noise = np.sqrt(noise_power / 2) * (np.random.normal(0, 1, dim) +
                                            1j * np.random.normal(0, 1, dim))
    else:
        noise = np.sqrt(noise_power) * np.random.normal(0, 1, dim)

    sig_noisy = sig + noise

    return sig_noisy

def razaoMatrizR(Y, L, M=128, SNR=20, potenciaSinal=1, retornarApenasParteReal=True, verbose=False):
    
    # OBTENDO A POTENCIA DO RUIDO N0
    N0 = potenciaSinal/(10**(SNR/10))
    
    # CRIANDO A MATRIZ IDENTIDADE E FAZENDO A OPERACAO HERMITIANA EM Y
    I  = np.identity(L)
    YH = np.conjugate(Y).T
    
    # CONSTRUINDO A MATRIZ R
    R = np.matmul(YH, Y)/M - N0 * I
    
    # DOIS MAIORES AUTOVALORES
    autovalores = np.sort(np.linalg.eigvals(R))
    a1 = autovalores[-1]
    a2 = autovalores[-2]
    
    if verbose:
        print("R:", np.matrix(R))
        print("\nAutovalores:", autovalores)
        print("\nMaior autovalor:", a1)
        print("\nSegundo maior autovalor:", a2)
        print("\nRazão:", (a1/a2).real)
    
    # RETORNO A RAZAO
    if retornarApenasParteReal:
        return (a1/a2).real
    else:
        return a1/a2