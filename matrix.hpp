#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <cassert>
#include <tuple>
#include <cmath>

using namespace std;
using std::size_t;

const double eps = 1e-8;

namespace sjtu {
    
    template <class T>
    class Matrix {
    private:
        size_t n, m;
        T *Mat;

    //construct & destroy
    public:
        Matrix() : n(0), m(0), Mat(nullptr) { }
        Matrix(size_t _n, size_t _m, T _init = T()) : n(_n), m(_m) {
            if (_n < 0 || _m < 0)
                throw invalid_argument("cannot be nagative");
            Mat = new T[n * m];
            for (size_t i = 0; i < n * m; i++)
                Mat[i] = _init;
        }
        explicit Matrix(std::pair<size_t, size_t> sz, T _init = T()) : Matrix(sz.first, sz.second, _init) { }

		Matrix(const Matrix &o) {
            n = o.rowLength(), m = o.columnLength();
            Mat = new T[n * m];

            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) = o(i, j);
		}
        template <class U>
        Matrix(const Matrix<U> &o) {
            n = o.rowLength(), m = o.columnLength();
            Mat = new T[n * m];
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) = o(i, j);
        }

        Matrix(Matrix &&o) noexcept : n(o.n), m(o.m), Mat(o.Mat) {
            o.Mat = nullptr;
        }

        Matrix(std::initializer_list<std::initializer_list<T>> il) {
            n = il.size(), m = il.begin()->size();
            Mat = new T[n * m];

            size_t p = 0, q;
            for (auto row : il) {
                if (row.size() != m) {
                    delete[] Mat;
                    throw invalid_argument("not a matrix");
                }
                q = 0;
                for (auto element : row)
                    (*this)(p, q++) = element;
                p++;
            }
        }

        ~Matrix() {
            delete[] Mat;
        }

        

    public:
        // assign matrix
        Matrix &operator=(const Matrix &o) {
            if (this == &o)
                return *this;
            delete[] Mat;

            n = o.n, m = o.m;
            Mat = new T[n * m];

            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) = o(i, j);
            return *this;
        }

        // assign
        template <class U>
        Matrix &operator=(const Matrix<U> &o) {
            //if (this == &o)
            //    return *this;
            delete[] Mat;

            n = o.rowLength(), m = o.columnLength();
            Mat = new T[n * m];

            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) = o(i, j);
            return *this;
        }
        

        Matrix &operator=(Matrix &&o) noexcept {
            if (this == &o)
                return *this;
            delete[] Mat;

            n = o.n, m = o.m, Mat = o.Mat;
            o.Mat = nullptr;
            return *this;
        }
        
        
        
    public:
        bool empty() const {
            return !n || !m;
        }
        size_t rowLength() const {
            return n;
        }
        size_t columnLength() const {
            return m;
        }
        std::pair<size_t, size_t> size() const {
            return std::make_pair(n, m);
        }

        void resize(size_t _n, size_t _m, T _init = T()) {
            if (_n == 0 || _m == 0) {
                delete[] Mat;
                n = 0, m = 0;
                return;
            }
            if (_n < 0 || _m < 0)
                throw invalid_argument("cannot be nagative");

            if (_n == n && _m == m)
                return;

            if (_n * _m == n * m) {
                n = _n, m = _m;
                return;
            }
            T *Mat1 = new T[_n * _m];

            for (size_t i = 0; i < _n * _m; i++)
                Mat1[i] = i < n * m ? Mat[i] : _init;
            delete[] Mat;
            Mat = Mat1;
            n = _n, m = _m;
        }
        void resize(std::pair<size_t, size_t> sz, T _init = T()) {
            resize(sz.first, sz.second, _init);
        }

        void clear() {
            delete[] Mat;
            Mat = nullptr;
            n = m = 0;
        }
        
    public:
        T &operator()(size_t i, size_t j) {
            if (i < 0 || i >= n || j < 0 || j >= m)
                throw std::invalid_argument("out of range");
            return Mat[i * m + j];
        }
        const T &operator()(size_t i, size_t j) const {
            if (i < 0 || i >= n || j < 0 || j >= m)
                throw std::invalid_argument("out of range");
            return Mat[i * m + j];
        }
        
        Matrix<T> row(size_t id) const {
            if (id < 0 || id >= n)
                throw std::invalid_argument("out of range");
            Matrix res(1, m);
            for (size_t i = 0; i < m; i++)
                res(0, i) = (*this)(id, i);
            return res;
        }
        Matrix<T> column(size_t id) const {
            if (id < 0 || id >= m)
                throw std::invalid_argument("out of range");
            Matrix res(n, 1);
            for (size_t i = 0; i < n; i++)
                res(i, 0) = (*this)(i, id);
            return res;
        }
        
        
    public:
        bool operator==(const Matrix &o) const {
            if (n != o.rowLength() || m != o.columnLength())
                return 0;
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    if (fabs((*this)(i, j) - o(i, j)) > eps)
                        return 0;
            return 1;
        }

        template <class U>
        bool operator==(const Matrix<U> &o) const {
            if (n != o.rowLength() || m != o.columnLength())
                return 0;
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    if (fabs((*this)(i, j) - o(i, j)) > eps)
                        return 0;
            return 1;
        }

        template <class U>
        bool operator!=(const Matrix<U> &o) const {
            return !(*this == o);
        }
        
        Matrix operator-() const {
            Matrix res = *this;
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    res(i, j) = -res(i, j);
            return res;
        }
        
        template <class U>
        Matrix &operator+=(const Matrix<U> &o) {
            if (n != o.n || m != o.m)
                throw std::invalid_argument("cannot add");
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) += o(i, j);
            return *this;
        }
        
        template <class U>
        Matrix &operator-=(const Matrix<U> &o) {
            if (n != o.n || m != o.m)
                throw std::invalid_argument("cannot subtract");
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) -= o(i, j);
            return *this;
        }
        
        template <class U>
        Matrix &operator*=(const U &x) {
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    (*this)(i, j) *= x;
            return *this;
        }
        
        Matrix tran() const {
            Matrix res(m, n);
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < m; j++)
                    res(j, i) = (*this)(i, j);
            return res;
        }
        
    public: // iterator
        class iterator {
		public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type        = T;
            using pointer           = T *;
            using reference         = T &;
            using size_type         = size_t;
            using difference_type   = std::ptrdiff_t;
            
            iterator() : p(nullptr), M(0), pos(0), U(0), L(0), D(0), R(0) {}
            iterator(const iterator &) = default;
            iterator &operator=(const iterator &) = default;
            iterator(const pointer &_p, const size_type &_M, const size_type &_pos, const size_t &_U, const size_t &_L,
                     const size_t &_D, const size_t &_R) : p(_p), M(_M), pos(_pos), U(_U), L(_L), D(_D), R(_R) { }
                
        private:
            pointer p;
            size_type M, pos;
            size_type U, L, D, R;

            std::pair<size_t, size_t> calc_pos() {
                return std::make_pair(pos / M, pos % M);
            }

            bool check_end_matrix(pointer p) {
                std::pair<size_t, size_t> pr = calc_pos();
                size_t x, y;
                std::tie(x, y) = pr;
                //x = pr.first, y = pr.second;
                return x == D && y == R;
            }

            bool check_in_matrix(pointer p) {
                std::pair<size_t, size_t> pr = calc_pos();
                size_t x, y;
                //std::tie(x, y) = pr;
                x = pr.first, y = pr.second;
                return x >= U && x <= D && y >= L && y <= R;
            }

        public:
            difference_type operator-(const iterator &o) {
                pointer s = p < o.p ? p : o.p;
                pointer t = p > o.p ? p : o.p;
                difference_type res = 0;
                while (s != t)
                    if (check_in_matrix(++s))
                        res++;
                return p > o.p ? res : -res;
            }
            
            iterator &operator+=(difference_type offset) {
                while (offset--)
                    p++, pos++;
                return *this;
            }
            
            iterator operator+(difference_type offset) const {
                iterator res = *this;
                res += offset;
                return res;
            }
            
            iterator &operator-=(difference_type offset) {
                while (offset--)
                    p--, pos--;
                return *this;
            }
            
            iterator operator-(difference_type offset) const {
                iterator res = *this;
                res -= offset;
                return res;
            }
            
            iterator &operator++() {
                do {
                    if (check_end_matrix(p)) {
                        p++, pos++;
                        break;
                    }
                    p++, pos++;
                } while (!check_in_matrix(p));
                return *this;
            }
            
            iterator operator++(int) {
                iterator res = *this;
                do {
                    if (check_end_matrix(p)) {
                        p++, pos++;
                        break;
                    }
                    p++, pos++;
                } while (!check_in_matrix(p));
                return res;
            }
            
            iterator &operator--() {
                do
                    p--, pos--;
                while (!check_in_matrix(p));
                return *this;
            }
            
            iterator operator--(int) {
                iterator res = *this;
                do
                    p--, pos--;
                while (!check_in_matrix(p));
                return *this;
            }
            
            reference operator*() const {
                return *p;
            }
            
            pointer operator->() const {
                return &(*p);
            }
            
            bool operator==(const iterator &o) const {
                return p == o.p;
            }
            
            bool operator!=(const iterator &o) const {
                return p != o.p;
            }
        };

        iterator begin() {
            return iterator(Mat, m, 0, 0, 0, n - 1, m - 1);
        }
        
        iterator end() {
            return iterator(Mat + n * m, m, n * m, 0, 0, n - 1, m - 1);
        }

        std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> L, std::pair<size_t, size_t> R) {
            size_t u, l, d, r;
            //std::tie(u, l) = L, std::tie(d, r) = R;
            u = L.first, l = L.second;
            d = R.first, r = R.second;
            if (u > d || l > r || u >= n || d < 0 || l >= m || r < 0)
                throw invalid_argument("out of range");
            iterator res1(Mat + u * m + l, m, u * m + l, u, l, d, r);
            iterator res2(Mat + d * m + r + 1, m, d * m + r + 1, u, l, d, r);
            return std::make_pair(res1, res2);
        }

    };
}

//
namespace sjtu
{
    template <class T, class U>
    auto operator*(const Matrix<T> &mat, const U &x) -> Matrix<decltype(mat(0, 0) + x)> {
        size_t n = mat.rowLength(), m = mat.columnLength();
        Matrix<decltype(mat(0, 0) + x)> res(n, m);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                res(i, j) = mat(i, j) * x;
        return res;
    }
    
    template <class T, class U>
    auto operator*(const U &x, const Matrix<T> &mat) -> Matrix<decltype(mat(0, 0) + x)> {
        return mat * x;
    }
    
    template <class U, class V>
    auto operator*(const Matrix<U> &a, const Matrix<V> &b) -> Matrix<decltype(a(0, 0) * b(0, 0))> {
        size_t n = a.rowLength(), k = a.columnLength(), m = b.columnLength();
        if (a.columnLength() != b.rowLength())
            throw std::invalid_argument("cannot multiply");
        Matrix<decltype(a(0, 0) * b(0, 0))> res(n, m);

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                for (size_t t = 0; t < k; t++)
                    res(i, j) += a(i, t) * b(t, j);
        return res;
    }

    template <class U, class V>
    auto operator+(const Matrix<U> &a, const Matrix<V> &b) -> Matrix<decltype(a(0, 0) + b(0, 0))> {
        size_t n = a.rowLength(), m = a.columnLength();
        if (n != b.rowLength() || m != b.columnLength())
            throw std::invalid_argument("cannot add");
        Matrix<decltype(a(0, 0) + b(0, 0))> res(n, m);

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                res(i, j) = a(i, j) + b(i, j);
        return res;
    }

    template <class U, class V>
    auto operator-(const Matrix<U> &a, const Matrix<V> &b) -> Matrix<decltype(a(0, 0) - b(0, 0))> {
        size_t n = a.rowLength(), m = a.columnLength();
        if (n != b.rowLength() || m != b.columnLength())
            throw std::invalid_argument("cannot subtract");
        Matrix<decltype(a(0, 0) - b(0, 0))> res(n, m);

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                res(i, j) = a(i, j) - b(i, j);
        return res;
    }    
}

#endif //SJTU_MATRIX_HPP
