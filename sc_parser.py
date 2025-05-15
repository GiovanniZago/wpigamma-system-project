from typing import List
from itertools import islice

import config
import puppi
import egamma

def generate_aiesim_file(puppis: List[dict], 
                         plio_width: int, 
                         egammas: List[dict] = None, 
                         pad_puppi: int = None, 
                         pad_egamma: int = None, 
                         file_out: str = None, 
                         debug: bool = False):
    
    """
    Check validity of plio_width
    """
    if plio_width not in [32, 64, 128]:
        raise ValueError("plio_width must be 32, 64 or 128")
    
    values_per_row = plio_width // 16 # suppose that each feature is 16-bit wide

    """
    Manage output file
    """
    fout = None
    if file_out:
        fout = open(config.AIE_DATA + f"/{file_out}.txt", "w")

    """
    Process puppi candidates
    """
    n_puppi = len(puppis)
    puppi_features = puppis[0].keys()

    if pad_puppi and (n_puppi < pad_puppi):
        for _ in range(pad_puppi - n_puppi):
            puppis.append(
                {key: 0 for key in puppi_features}
            ) 

        n_puppi = pad_puppi

    elif n_puppi > pad_puppi:
        raise ValueError("n_puppi is bigger than pad_puppi")

    """ 
    Process egamma candidates, if any
    """
    if egammas:
        n_egamma = len(egammas)
        egamma_features = egammas[0].keys()

        if pad_egamma and (n_egamma < pad_egamma):
            for _ in range(pad_egamma - n_egamma):
                egammas.append(
                    {key: 0 for key in egamma_features}
                )

            n_egamma = pad_egamma

        elif n_egamma > pad_egamma:
            raise ValueError("n_egamma is bigger than pad_egamma")
    
    """ 
    Flatten puppis and egammas
    """
    # remember that in multiple list comprehension like the following one, 
    # the outer "for" is the one changing faster (for can in puppis, in this case)
    puppi_values = [cand[key] for key in puppi_features for cand in puppis]
    egamma_values = [cand[key] for key in egamma_features for cand in egammas]
    
    values = puppi_values + egamma_values

    if len(values) % values_per_row:
        print("The total number of data values is not a multiple of the values_per_row")    
        return 1
    
    itervalues = iter(values)

    while True:
        batch = list(islice(itervalues, values_per_row))

        if not batch:
            break
        
        string = " ".join([str(val) for val in batch])

        if debug:
            print(string)

        if fout:
            fout.write(string + "\n")

    if fout:
        fout.close()


if __name__ == "__main__":
    # read puppis
    fbin = open(config.DATA + "/puppi_WPiGamma_PU200.dump", "rb")

    header_byte = 496
    _, _, puppis = puppi.get_puppi_cands(header_byte, fbin, ["eta", "pid"], debug=True)

    fbin.close()

    # read egammas
    fbin = open(config.DATA + "/egamma_WPiGamma_PU200.dump", "rb")

    header_byte = 344
    _, _, egammas = egamma.get_egamma_cands(header_byte, fbin, ["eta"], debug=True)

    fbin.close()

    # get aiesim file
    generate_aiesim_file(puppis, 32, egammas, 32, 32, file_out="foo1")
