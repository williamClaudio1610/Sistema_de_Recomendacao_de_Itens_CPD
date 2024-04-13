	gcc -o projectoSerial.o projectoSerial.c -c -g -fopenmp
	gcc -o projectoSerial *.o -fopenmp
 run:
	./projectoSerial

clean:
	rm *.o
	rm projectoSerial
