# Multiframe Motion Coupling for Video Super Resolution


## Authors
* Hendrik Dirks ([hendrik.dirks@wwu.de](mailto:hendrik.dirks@wwu.de))*
* Jonas Geiping ([jonas.geiping@uni-siegen.de](jonas.geiping@uni-siegen.de))*
* Michael Moeller ([michael.moeller@uni-siegen.de](michael.moeller@uni-siegen.de))*

### Dependencies
* MATLAB 
* CUDA

## Please make sure to install all submodules and check their installation instructions

### Quick start
for a quick installation on a Linux system run the following to install all submodules
```
git clone https://github.com/HendrikMuenster/superResolution
cd superResolution
git submodule update --init
cd flexBox
git submodule update --init # we only need to initialize this here
mkdir flexBox_CPP/build
cd flexBox_CPP/build
cmake -DUSE_CUDA=ON -DUSE_OPENMP=ON ../source
make && make install
cd ../../../prost
mkdir prost/build
cd prost/build
cmake ..
make
mkdir build
cd prost/build
cmake ..
make
```
If necessary, run cmake-gui instead of cmake to easily set manual options (e.g. MATLAB installation path).
It is also possible to install flexBox without CUDA (or without compiling the C++ submodule at all and
running just MATLAB code. Please notify us, if you have trouble compiling one of these modules.


## Usage
For a quick example, given a video as a 4-D matlab array 'videoLowRes' that we want to upsample by a factor of 4,
run the following commands in MATLAB:
```
MMC = MultiframeMotionCoupling(videoLowRes,'zoom_factor',4);
MMC.init;
MMC.run;
videoHighRes = MMC.result1;
```

Check the testfile.m for further instructions.

## Options

To reproduce the results of the paper, use the (standard) parameters
```
MMC.framework = 'prost';
MMC.comp_mode = 'accurate';
```
However you can also optionally use the 'flexBox' framework for the super resolution step. 
Other than that, can change to computation mode to 'faster' (reducing speed) or 'fastest' 
(reducing memory consumption and increasing speed), although no guarantee can be given on the
quality of the results.

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