#ifndef makevolume_core_h
#define makevolume_core_h

// functor for indexing
template<class T>
struct index_sorting {
  const T _values;

public:
  index_sorting(const T values) : _values(values) {};
  bool operator() (long unsigned a, long unsigned b) const {
    return _values[a] < _values[b]; }    
};

long unsigned* indexing(long unsigned n, float* array);

#endif // makevolume_core_h
