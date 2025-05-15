import sys
import numpy as np
import config

HEADER_SIZE = 8
WORD_SIZE   = 8

def unpack_puppi(byte_line, reals=False):
    data = int.from_bytes(byte_line, sys.byteorder)

    pid         = (data >> 37) & 0b111          # Bits 39-37 (3 bits)
    phi_sign    = (data >> 36) & 0b1            # Bit 36 is the sign bit for phi
    phi_payload = (data >> 26) & 0x3FF          # Bits 35-26 (10 bits) are the non-sign bits for phi
    eta_sign    = (data >> 25) & 0b1            # Bit 25 is the sign bit for eta
    eta_payload = (data >> 14) & 0x7FF          # Bits 24-14 (11 bits) are the non-sign bits for eta
    pt          = data & 0x3FFF                 # Bits 13-00 (14 bits) 

    phi = 0
    for ii in range(phi_payload.bit_length()):
        cur_bit = (phi_payload >> ii) & 0b1

        if cur_bit:
            phi += 2 ** ii

    phi += (-1) * phi_sign * 2 ** (11 - 1)

    eta = 0
    for ii in range(eta_payload.bit_length()):
        cur_bit = (eta_payload >> ii) & 0b1

        if cur_bit:
            eta += 2 ** ii

    eta += (-1) * eta_sign * 2 ** (12 - 1)

    if reals:
        return {
            "pt": pt * 0.25,
            "eta": eta * (np.pi / 720), 
            "phi": phi * (np.pi / 720), 
            "pid": pid 
        }
    
    else:
        return {
            "pt": pt,
            "eta": eta, 
            "phi": phi, 
            "pid": pid 
        }
    
def get_puppi_cands(header_byte, file, features=None, reals=False, debug=False):
    file.seek(header_byte)

    if features is None:
        features = ["pt", "eta", "phi", "pid"]

    foo = int.from_bytes(file.read(HEADER_SIZE), sys.byteorder)

    vld_header = (foo >> 62) & 0b11                    # Bits 63-62 (2 bits)
    err_bit    = (foo >> 61) & 0b1                     # Bit 61 (1 bit)
    lr_num     = (foo >> 56) & 0b11111                 # Bits 60-56 (5 bits)
    orbit_cnt  = (foo >> 24) & 0xFFFFFFFF              # Bits 55-24 (32 bits)
    bx_cnt     = (foo >> 12) & 0xFFF                   # Bits 23-12 (12 bits)
    n_cands    = foo & 0xFF                            # Bits 11-00 (12 bits) BUT ONLY 8 EFFECTIVELY USED (look at the binary mask indeed)

    next_header_byte = header_byte + (n_cands + 1) * WORD_SIZE

    print(f"Puppi header byte: {header_byte}, vld_header: {vld_header}, err_bit: {err_bit}, "
            f"lr_num: {lr_num}, orbit_cnt: {orbit_cnt}, bx_cnt: {bx_cnt}, n_cand: {n_cands}, "
            f"next_header_byte: {next_header_byte}")
    
    data = []

    if vld_header != 2:
        print("Corrupted event: vld_header != 0b10")
        return 1, next_header_byte, data

    if err_bit == 1:
        print("Corrupted event: err_bit = 0b1")
        return 1, next_header_byte, data
    
    for _ in range(n_cands):
        particle = unpack_puppi(file.read(WORD_SIZE), reals)
        data.append(
            {key: particle[key] for key in features}
        )
    
    if debug:
        for particle in data:
            print(particle)

    return 0, next_header_byte, data
    
if __name__ == "__main__":
    file = open(config.DATA + "/puppi_WPiGamma_PU200.dump", "rb")
    get_puppi_cands(496, file, reals=False, debug=True)