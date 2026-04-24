# code-for-mobility-spectrum-method

code for mobility spectrum method


The related paper is:

1. Ting-Na Shao, Zi-Tao Zhang, Yu-Jie Qiao, Qiang Zhao, Hai-Wen Liu*, Xin-xiang Chen, Wei-Min Jiang, Chun-Li Yao, Xing-Yu Chen, Mei-Hui Chen, Rui-Fen Dou, Chang-Min Xiong, Guang-Ming Zhang*, Yi-Feng Yang*, Jia-Cai Nie*, Kondo scattering in underdoped Nd1−xSrxNiO2 infinite-layer superconducting thin films, National Science Review, 10(11), nwad112 (2023), DOI: https://doi.org/10.1093/nsr/nwad112.

2. Jianchao Meng, Xinxiang Chen, Tingna Shao, Mingrui Liu, Weimin Jiang, Zitao Zhang, Changmin Xiong*, Ruifen Dou*, Jiacai Nie*, Doping-enhanced robustness of anomaly-related magnetoresistance in WTe2±α flakes, Chinese Physics B, 32(4), 047502 (2023), DOI: https://doi.org/10.1088/1674-1056/acb423.

3. Jianchao Meng, Xinxiang Chen, Mingrui Liu, Weimin Jiang, Zhe Zhang, Jingzhuo Ling, Tingna Shao, Chunli Yao, Lin He*, Ruifen Dou*, Changmin Xiong* and Jiacai Nie*, Large linear magnetoresistance caused by disorder in WTe2–δ thin film, Journal of Physics: Condensed Matter 32, 355703 (2020), DOI: https://doi.org/10.1088/1361-648X/ab8d74.


# Mobility Spectrum Analysis

This program reconstructs the mobility distribution \(p(\mu)\) from magnetic-field-dependent conductivity data using a maximum-entropy-like iterative method.

The input data are longitudinal and Hall conductivity components measured at different magnetic fields. The program outputs the reconstructed mobility spectrum and estimates the dominant negative- and positive-mobility carrier contributions.



# How to use it

Some modules need to involve:

need FFTW

For example:

need \sigma_0.txt and b_segma_100K.txt


gcc main.c -o mobility -lfftw3 -lm
./mobility 100

