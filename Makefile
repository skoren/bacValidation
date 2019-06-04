main: bamcat.o samToErrorRate.o
	g++ -o samToErrorRate -lz -lpthread bamcat.o samToErrorRate.o
	rm -f bamcat.o
	rm -f samToErrorRate.o

bamcat.o: bamcat.c
	gcc -g -c bamcat.c

samToErrorRate.o: samToErrorRate.C
	g++ -g -c samToErrorRate.C

clean:
	rm -f bamcat.o
	rm -f samToErrorRate.o
	rm -f samToErrorRate
