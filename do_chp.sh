cc -O4  -I/Users/bschill/rtl/rtl-sdr/include -I/opt/local/include/libusb-1.0  -c    chp.c -o chp.c.o
cc  -O3 -DNDEBUG -Wl,-search_paths_first -Wl,-headerpad_max_install_names   chp.c.o  -o rtl_chp  ../build/src/librtlsdr.0.5git.dylib /opt/local/lib/libusb-1.0.dylib -lpthread -lm 

