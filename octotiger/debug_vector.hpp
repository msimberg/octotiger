/*
 * debug_vector.hpp
 *
 *  Created on: Jan 30, 2019
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_DEBUG_VECTOR_HPP_
#define OCTOTIGER_DEBUG_VECTOR_HPP_

#include <array>
#include <vector>
#include <typeinfo>
#include <iostream>


template<class T, int N>
struct debug_array: public std::array<T,N> {
	using A = std::array<T,N>;
	using value_type = typename A::value_type;
	using A::A;
	const auto& operator[](int i) const {
		check_range(i);
		return A::operator[](i);
	}

	auto& operator[](int i) {
		check_range(i);
		return A::operator[](i);
	}

	explicit debug_array(const std::array<T,N>& other) :
			A(other) {
	}

	debug_array(std::initializer_list<T> l ) {
		std::copy(l.begin(),l.end(), A::begin());
	}



private:
	void check_range(int i) const {
		constexpr auto lb = 0;
		const auto ub = N;
		if (i < lb || i >= ub) {
			std::cout << "INDEX OUT OF RANGE type = ";
			std::cout << typeid(debug_array).name();
			std::cout << " Index is " << std::to_string(i);
			std::cout << " it should be between ";
			std::cout << std::to_string(lb) << " and " << std::to_string(ub - 1) << ".\n";
			std::abort();
		}
	}

	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & *dynamic_cast<A*>(this);
	}
};

template<class T>
struct debug_vector: public std::vector<T> {
	using A = std::vector<T>;
	using value_type = typename A::value_type;
	using A::A;
	const auto& operator[](int i) const {
		check_range(i);
		return A::operator[](i);
	}

	auto& operator[](int i) {
		check_range(i);
		return A::operator[](i);
	}

	explicit debug_vector(const std::vector<T>& other) :
			A(other) {
	}

	explicit debug_vector(std::vector<T>&& other) :
			A(std::move(other)) {
	}

private:
	void check_range(int i) const {
		constexpr auto lb = 0;
		const auto ub = A::size();
		if (i < lb || i >= ub) {
			std::cout << "INDEX OUT OF RANGE type = ";
			std::cout << typeid(debug_vector).name();
			std::cout << " Index is " << std::to_string(i);
			std::cout << " it should be between ";
			std::cout << std::to_string(lb) << " and " << std::to_string(ub - 1) << ".\n";
			std::abort();
		}
	}

	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & *dynamic_cast<A*>(this);
	}
};

template<>
struct debug_vector<bool> : public std::vector<bool> {
	using A = std::vector<bool>;
	using value_type = typename A::value_type;
	using A::A;
	auto operator[](int i) const {
		check_range(i);
		return A::operator[](i);
	}

	auto operator[](int i) {
		check_range(i);
		return A::operator[](i);
	}
	debug_vector(const std::vector<bool>& other) :
			A(other) {
	}

	debug_vector(std::vector<bool>&& other) :
			A(std::move(other)) {
	}


	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & *dynamic_cast<A*>(this);
	}

private:
	void check_range(int i) const {
		constexpr auto lb = 0;
		const auto ub = A::size();
		if (i < lb || i >= ub) {
			std::cout << "INDEX OUT OF RANGE type = ";
			std::cout << typeid(debug_vector).name();
			std::cout << " Index is " << std::to_string(i);
			std::cout << " it should be between ";
			std::cout << std::to_string(lb) << " and " << std::to_string(ub - 1) << ".\n";
			std::abort();
		}
	}
};
namespace oct {
#ifdef NDEBUG
template<class T>
using vector = std::vector<T>;

template<class T, int N>
using array = std::array<T,N>;

#else
template<class T>
using vector = debug_vector<T>;

template<class T, int N>
using array = debug_array<T,N>;
#endif
}
#endif /* OCTOTIGER_DEBUG_VECTOR_HPP_ */
