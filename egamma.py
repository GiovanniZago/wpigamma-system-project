import sys
import config 
import numpy as np

HEADER_SIZE = 8
LOW_SIZE    = 8
HIGH_SIZE   = 4

def unpack_egamma(byte_line, reals=False):
    """unpack_egamma

    Receives an egamma byte string and parses the egamma features 
    into a dictionary.

    Parameters
    ----------
    byte_line : bytes
        12-byte word containing the features of an egamma candidate
    reals : bool, optional
        If true, the egamma features are returned in the dictionary as 
        floating point values instead of integers, by default False

    Returns
    -------
    dict
        Dictionary containing parsed egamma features, which are
        pt, eta, phi and quality.
    """
    data = int.from_bytes(byte_line, sys.byteorder)

    # pt: bits 1–16, unsigned
    pt = (data >> 1) & 0xFFFF

    # phi: bits 17–29, sign bit at 29
    phi_sign = (data >> 29) & 0b1
    phi_payload = (data >> 17) & 0xFFF

    # eta: bits 30–43, sign bit at 43
    eta_sign = (data >> 43) & 0b1
    eta_payload = (data >> 30) & 0x1FFF

    # quality: bits 44–47, unsigned
    quality = (data >> 44) & 0xF

    phi = 0
    for ii in range(phi_payload.bit_length()):
        cur_bit = (phi_payload >> ii) & 0b1

        if cur_bit:
            phi += 2 ** ii

    phi += (-1) * phi_sign * 2 ** (13 - 1)


    eta = 0
    for ii in range(eta_payload.bit_length()):
        cur_bit = (eta_payload >> ii) & 0b1

        if cur_bit:
            eta += 2 ** ii

    eta += (-1) * eta_sign * 2 ** (14 - 1)

    if reals:
        return {
            "pt": pt * 0.03125,
            "phi": phi * np.pi / (2 ** 12),
            "eta": eta * np.pi / (2 ** 12), 
            "quality": quality
        }
    
    else:
        return {
            "pt": pt,
            "phi": phi,
            "eta": eta, 
            "quality": quality
        }

def get_egamma_event(header_byte, file, file_out=None, features=None, pad_to=None, reals=False, debug=False):
    """get_egamma_event

    Allows for inspecting the data of a binary file

    Parameters
    ----------
    header_byte : int
        Index of the byte corresponding to the header of an event
    file : open file object
        File from which the data has to be read
    file_out : str, optional
        Name of the file in which the unpacked data is written to, by default None
    features : list, optional
        List of the features to be print to file, if present
    pad_to : int, optional
        Number of desided candidates in the current event. If n_cands is less than
        pad_to, then the due number of trailing zeros is added to each feature.
    reals : bool, optional
        Convert unpacked data to floating point, by default False
    debug : bool, optional
        Print particle data, by default False

    Returns
    -------
    is_corrupt : int
        1 if the data in the input file is corrupted, otherwise 0
    next_header_byte : int
        Index of the byte at which the next event header is located
    """
    file.seek(header_byte)

    if features is None:
        features = ["pt", "eta", "phi", "quality"]

    """
    Read header
    """
    foo = int.from_bytes(file.read(HEADER_SIZE), sys.byteorder)
    
    """ 
    Unpack header. Pay attention that the words intended in
    n_words are 64-bit (8-byte) words, as it is the atomic
    byte word of the data format. A candidate is split between
    two words (64-bit + 32-bit) so, since each candidate must 
    be complete, the minimum unit that one has to parse is 
    three blocks of 64-bit words. 
    """
    n_words          = foo & 0xFFF
    n_cands          = (n_words * 2) // 3
    next_header_byte = header_byte + (n_words + 1) * 8

    print(f"Header byte: {header_byte}, n_words: {n_words}, n_cands: {n_cands}, next_header_byte: {next_header_byte}")

    # calculate 12-byte blocks that hav
    inverse          = False
    n_blocks         = n_words // 3
    spare_words      = n_words % 3

    if (spare_words == 1) or (n_words == 0):
        return 1, next_header_byte

    data = []

    for _ in range(n_blocks):
        for _ in range(2):
            if inverse:
                bytes_high = file.read(HIGH_SIZE)
                bytes_low = file.read(LOW_SIZE)
            else:
                bytes_low = file.read(LOW_SIZE)
                bytes_high = file.read(HIGH_SIZE)

            bytes_tot = bytes_low + bytes_high
            particle = unpack_egamma(bytes_tot, reals)

            data.append(
                {key: particle[key] for key in features}
            )

            inverse = not inverse

    if spare_words:
        if inverse:
            bytes_high = file.read(HIGH_SIZE)
            bytes_low = file.read(LOW_SIZE)
        else:
            bytes_low = file.read(LOW_SIZE)
            bytes_high = file.read(HIGH_SIZE)

        bytes_tot = bytes_low + bytes_high
        particle = unpack_egamma(bytes_tot, reals)

        data.append(
            {key: particle[key] for key in features}
        )

    if debug:
        for particle in data:
            print(particle)

    if file_out:
        fout = open(config.AIE_DATA + f"/{file_out}.txt", "w")

        for particle in data:
            fout.write(" ".join([str(particle[key]) for key in particle]) + "\n")

        if pad_to and (n_cands < pad_to):
            for _ in range(pad_to - n_cands):
                fout.write(" ". join(["0" for _ in features]) + "\n")

        fout.close()
    
    return 0, next_header_byte

if __name__ == "__main__":
    fbin = open(config.DATA + "/egamma_WPiGamma_PU200.dump", "rb")

    header_byte = 0
    get_egamma_event(header_byte, fbin, "in1", ["phi", "quality"], pad_to=32)


    fbin.close()