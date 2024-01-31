/*
 * File: ./Extern/compare_vector.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#ifndef COMPARE_VECTOR_H
#define COMPARE_VECTOR_H
#include <iostream>
#include <vector>

template <class T>
void compare_num_data(const std::vector<T> vector1,
    const std::vector<T> vector2)
{
    if (vector1.size() < vector2.size()) {
        std::cout << "\""
                  << "vector1"
                  << "\" has only " << vector1.size() << " elements while \""
                  << "vector2"
                  << "\" has " << vector2.size()
                  << ". Can only compare vectors of equal size." << std::endl;
    } else if (vector2.size() < vector1.size()) {
        std::cout << "\""
                  << "vector2"
                  << "\" has only " << vector2.size() << " elements while \""
                  << "vector1"
                  << "\" has " << vector1.size()
                  << ". Can only compare vectors of equal size." << std::endl;
    } else {
        unsigned vector1_smaller_count = 0;
        unsigned vector2_smaller_count = 0;

        unsigned first_difference = (unsigned)-1;

        for (unsigned i = 0; i < vector1.size(); ++i) {
            if (vector1[i] < vector2[i])
                ++vector1_smaller_count;
            if (vector2[i] < vector1[i])
                ++vector2_smaller_count;
            if (vector1[i] != vector2[i] && first_difference == (unsigned)-1)
                first_difference = i;
        }

        if (vector1_smaller_count == 0 && vector2_smaller_count == 0) {
            std::cout << "The vectors are the same and have " << vector1.size()
                      << " elements." << std::endl;
        } else {
            std::cout << "The vectors differ. " << vector1_smaller_count
                      << " elements are smaller in \""
                      << "vector1"
                      << "\". " << vector2_smaller_count
                      << " elements are smaller in \""
                      << "vector2"
                      << "\". "
                      << (vector1.size() - vector1_smaller_count - vector2_smaller_count)
                      << " elements are the same. "
                      << (vector1_smaller_count + vector2_smaller_count)
                      << " elements are different. "
                      << "The vectors have " << vector1.size() << " elements."
                      << "The first element that differ is at index "
                      << first_difference << "." << std::endl;
        }
    }
}

#endif