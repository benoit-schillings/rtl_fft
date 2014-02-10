cc -O4  -I/Users/bschill/rtl/rtl-sdr/include -I/opt/local/include/libusb-1.0  -c    fft.c -o fft.c.o
cc  -O3 -DNDEBUG -Wl,-search_paths_first -Wl,-headerpad_max_install_names   fft.c.o  -o rtl_fft  ../build/src/librtlsdr.0.5git.dylib /opt/local/lib/libusb-1.0.dylib -lpthread -lm 

