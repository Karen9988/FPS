# FPS_P/N
**FPS_P/N** is a two-dimensional mass spectrometry utilization program with precursor ion determination for accurate annotation of anthocyanins and flavonol glycosides.

### Prerequisite
#### 1. Msdial tool
You should fetch `msdial_tool` from [here](https://github.com/Karen9988/MsdialWorkbench/releases/) which is submodule of our program.

`msdial_tool` releases run on Windows with `.Net` Framework.

#### 2. python requirements
This Python program runs ok when all requirements installed by command:
`pip install -r requirements.txt`

### Usage
```shell
python fps_p_n.py --msdial_tool_path /path/to/msdial_tool \
                  --deconv_ms2_input_dir /path/to/ms2_mzML_dir \
                  --deconv_ms2_output_dir /path/to/ms2_deconvolution_result_dir \
                  --deconv_fps_ms1_input_dir /path/to/ms1_mzML_dir  \
                  --deconv_fps_ms1_output_dir /path/to/ms1_deconvolution_result_dir \
                  --networking_output_dir /path/to/molecularnetworking/result_dir \
                  --cluster_output_dir /path/to/cluster/output_dir  \
                  --msp_path /path/to/library/dir \
                  --annotation_output_dir  /path/to/annotation/output_dir
```
