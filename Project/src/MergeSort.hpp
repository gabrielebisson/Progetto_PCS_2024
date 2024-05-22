//--------------- MERGE SORT DECRESCENTE--------------------
#pragma once
#include <vector>
#include <algorithm> // per std::copy

namespace SortLibrary {

template<typename T>
void Merge(std::vector<T>& v,
           const unsigned int& sx,
           const unsigned int& cx,
           const unsigned int& dx){

    unsigned int i = sx;
    unsigned int j = cx + 1;

    std::vector<T> b;
    b.reserve(dx - sx + 1);

    while (i <= cx && j <= dx) {
        if (std::get<2>(v[i]) >= std::get<2>(v[j])) // Ordinamento decrescente in base alla lunghezza (terzo elemento della tupla)
            b.push_back(v[i++]);
        else
            b.push_back(v[j++]);
    }

    while (i <= cx)
        b.push_back(v[i++]);
    while (j <= dx)
        b.push_back(v[j++]);

    std::copy(b.begin(), b.end(), v.begin() + sx); // Copia i valori ordinati nel vettore originale
}

template<typename T>
void MergeSort(std::vector<T>& v, const unsigned int& sx, const unsigned int& dx) {
    if (sx < dx) {
        unsigned int cx = (sx + dx) / 2;
        MergeSort(v, sx, cx);
        MergeSort(v, cx + 1, dx);
        Merge(v, sx, cx, dx);
    }
}

template<typename T>
void MergeSort(std::vector<T>& v) {
    if (!v.empty()) {
        MergeSort(v, 0, v.size() - 1);
    }
}

}


