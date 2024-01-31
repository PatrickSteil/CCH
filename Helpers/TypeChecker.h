/*
 * File: ./Helpers/TypeChecker.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

template <typename T>
struct is_numerical {
    static constexpr bool value = std::is_integral<T>::value || std::is_floating_point<T>::value;
};