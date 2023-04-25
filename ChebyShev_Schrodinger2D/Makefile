FFTW_DIR=/opt/homebrew/

BC: SP2D_AdaptiveMultiChebySpectral_BC1.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP2D_AdaptiveMultiChebySpectral_BC1.cpp -o SP2D_AdaptiveMultiChebySpectral_BC1 --std=c++17  -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g -O3

AdaptiveMultiChebySpectral_exp: SP2D_AdaptiveMultiChebySpectral_exp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP2D_AdaptiveMultiChebySpectral_exp.cpp -o SP2D_AdaptiveMultiChebySpectral_expADI --std=c++17  -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g

exp1Dinterp: SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp -o SP2D_AdaptiveMultiChebySpectral_expADI1Dinterp --std=c++17  -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g

exp1Dinterp: SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp -o SP2D_AdaptiveMultiChebySpectral_expADI1Dinterp --std=c++17  -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g

exp1Dinterp2: SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp headers_Cheby/input.h headers_Cheby/write.h
	g++ SP2D_AdaptiveMultiChebySpectral_exp1Dinterp.cpp -o SP2D_AdaptiveMultiChebySpectral_expADI1Dinterp2 --std=c++17  -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -g
