CC = gcc
CFLAGS =-O3
CFLAGS1=-O3 -g -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99

all :
	$(CC) $(CFLAGS) -c LizarMong_KEM_Comfort.c main.c randombytes.c fips202.c xe5.c
	$(CC) $(CFLAGS) -o LizarMong_KEM_Comfort LizarMong_KEM_Comfort.o main.o randombytes.o fips202.o xe5.o
	
run : all
	./LizarMong_KEM_Comfort

clean :
	rm -f *.o
	rm -f LizarMong_KEM_Comfort

new :
	make clean
	make all
	./LizarMong_KEM_Comfort
