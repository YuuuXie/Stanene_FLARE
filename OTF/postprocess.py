import numpy as np
from flare.otf_parser import OtfAnalysis

filename = '../OTF_train/otf_run_100ps.out'
otf = OtfAnalysis(filename)
msds = otf.msds
dft_frames = np.array(otf.dft_frames)
np.save('MSD', msds)
np.save('DFT_calls', dft_frames)

