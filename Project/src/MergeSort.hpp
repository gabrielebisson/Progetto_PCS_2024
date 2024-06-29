//MergeSort decrescente leggermente modificato per le esigenze del progetto

#pragma once
#include <vector>
#include <algorithm> // per std::copy

namespace SortLibrary {

template<typename T1, typename T2>
void Merge(std::vector<T1>& v,
           const std::vector<T2>& rif,
           const unsigned int& sx,
           const unsigned int& cx,
           const unsigned int& dx){

    unsigned int i = sx;
    unsigned int j = cx + 1;

    std::vector<T1> b;
    b.reserve(dx - sx + 1);

    while (i <= cx && j <= dx) {
        if (rif[v[i]] >= rif[v[j]]) // Ordinamento decrescente in base alla lunghezza (terzo elemento della tupla)
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

template<typename T1, typename T2>
void MergeSort(std::vector<T1>& v,
               const std::vector<T2>& rif,
               const unsigned int& sx,
               const unsigned int& dx) {
    if (sx < dx) {
        unsigned int cx = (sx + dx) / 2;
        MergeSort(v,rif, sx, cx);
        MergeSort(v,rif, cx + 1, dx);
        Merge(v,rif, sx, cx, dx);
    }
}

template<typename T1, typename T2>
void MergeSort(std::vector<T1>& v,
               const std::vector<T2>& rif) {
    if (!v.empty()) {
        MergeSort(v,rif, 0, v.size() - 1);
    }
}

}


