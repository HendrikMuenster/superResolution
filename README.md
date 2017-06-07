# superResolution


## Multiframe Motion Coupling for Video Super Resolution

[This is preliminary code]

## Authors
* Hendrik Dirks ([hendrik.dirks@wwu.de](mailto:hendrik.dirks@wwu.de))*
* Jonas Geiping ([jonas.geiping@uni-siegen.de](jonas.geiping@uni-siegen.de))*
* Michael Moeller ([michael.moeller@uni-siegen.de](michael.moeller@uni-siegen.de))*

### Dependencies
* MATLAB 
* CUDA 

## Please make sure to install all submodules and check their installation instructions

### Quick start
for a quick installation on a Linux system run
```
git clone https://github.com/HendrikMuenster/superResolution
cd superResolution
git submodule --init --recursive
mkdir prost/build
cd prost/build
cmake ..
make -j8
cd ../../flexBox/flexBox_CPP/source
mkdir build
cd build 
cmake -DUSE_CUDA=ON -DUSE_OPENMP=ON ../
make -j8
make install
```
If necessary, run cmake-gui instead of cmake to easily set manual options (e.g. MATLAB installation path)


## Usage
For a quick example, given a video as a 4-D matlab array 'videoLowRes' that we want to upsample by a factor of 4,
run the following commands in MATLAB:
```
MMC = MultiframeMotionCoupling(videoLowRes);
MMC.factor = 4;
MMC.init;
MMC.run;
videoHighRes = MMC.result1;
```

Check the testfile.m for further instructions.

## Reference
While the article is under review, refer to
```
@Article{MMC,
  Title         = {Multiframe Motion Coupling via Infimal Convolution Regularization
                   for Video Super Resolution},
  Author        = { Hendrik Dirks Jonas Geiping and Daniel Cremers and Michael Moeller},
  Journal       = {ArXiv e-prints},
  url           = {http://arxiv.org/abs/1611.07767}
  Year          = {2016},
  Month         = dec,
  Primaryclass  = {cs.CV}
}
```