# libadf.a: kernels.h kernels.cpp graph.h graph.cpp


setup:
	mkdir aie hls
	source ~/xilinx_setup.sh
	xbmgmt examine


clean:
	rm -r aie hls