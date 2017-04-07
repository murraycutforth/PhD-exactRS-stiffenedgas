TARGET := main.exe
CFLAGS := -Wall -c -std=c++11 -O3 -g

$(TARGET) : main.o exact_RS_idealgas.o
	g++ $^ -O3 -o $(TARGET)

%.o: %.cpp
	g++ $(CFLAGS) -o $@ $<

clean:
	rm *.o $(TARGET)
