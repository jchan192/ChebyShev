FFTW_DIR=/opt/homebrew/

AdaptiveMultiChebySpectral: SP3D_AdaptiveMultiChebySpectral.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP3D_AdaptiveMultiChebySpectral.cpp -o SP3D_AdaptiveMultiChebySpectral --std=c++17 -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -O3

AdaptiveMultiChebySpectral_exp: SP3D_AdaptiveMultiChebySpectral_exp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP3D_AdaptiveMultiChebySpectral_exp.cpp -o SP3D_AdaptiveMultiChebySpectral_exp --std=c++17 -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g -O3

PGaussian_exp: SP3D_AdaptiveMultiChebySpectral_exp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP3D_AdaptiveMultiChebySpectral_exp.cpp -o SP3D_PGaussian --std=c++17 -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g -O3
