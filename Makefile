CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wshadow -Werror -fopenmp -O3 -DNDEBUG  -ffast-math -march=native
EXTRA = -ftree-vectorizer-verbose=1
TARGET =mgsolve
HXX= MG.h Timer.h

OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) $(HXX) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET)
