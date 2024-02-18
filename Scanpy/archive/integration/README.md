# Tools for data integration / batch correction

## Test time usage as a function of the number of libraries integrated

* [BBKNN](test_time_bbknn.py)
* [Harmony py](test_time_harmony.py)
* [scVI]()

## Test quality of integration

* [BBKNN]()
* [Harmony py]()
* [scVI]()


## Other

* Figure out [Geneformer]()

## Installation / resources

### scvi

* Tool available [here](https://docs.scvi-tools.org/en/stable/index.html)
* qrsh -l h="*epsilon*",nvgpu=4,mem_requested=80G,h_vmem=80G,tmp_requested=100G -now no
* micromamba create -n scvi-env_gpu python=3.9
* pip install torch torchvision torchaudio -f https://download.pytorch.org/whl/cu114/torch_stable.html
* micromamba install jax jaxlib -c conda-forge
* micromamba install scvi-tools -c conda-forge
* pip install chex==0.1.8
* pip install torch torchvision torchaudio -f https://download.pytorch.org/whl/cu114/torch_stable.html
* micromamba activate /home/anncuo/y/envs/scvi_gpu

### Harmony

* Harmony (original R version): https://portals.broadinstitute.org/harmony/articles/quickstart.html
* Harmony python implementation used here: https://github.com/slowkow/harmonypy ([example usage](https://support.parsebiosciences.com/hc/en-us/articles/7704577188500-How-to-analyze-a-1-million-cell-data-set-using-Scanpy-and-Harmony))

### BBKNN

* BBKNN: https://github.com/Teichlab/bbknn ([example usage](https://nbviewer.org/github/Teichlab/bbknn/blob/master/examples/simulation.ipynb))
