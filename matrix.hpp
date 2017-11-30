#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <exception>
#include <vector>
#include <iostream>
#include <utility>
#include <functional>
#include <initializer_list>
#include <iterator>

using std::size_t;


// matrix
namespace sjtu
{
	
	template <class T>
	class Matrix
	{
	private:
		T *storage = nullptr;
		size_t n = 0, m = 0;
		
	public:
		Matrix() = default;
		
		Matrix(size_t n, size_t m, T _init = T())
		{
			init(n, m, _init);
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
			init(sz.first, sz.second, _init);
		}
		
		// vector
		// The parameter _init is not provided, otherwise this ctor will be
		// confused with the ctor Matrix(size_t n, size_t m, T _init = T())
		explicit Matrix(size_t n)
		{
			init(1, n, T());
		}
		
		Matrix(const Matrix &o)
		{
			init(o.rowLength(), o.columnLength());
			size_t cnt = 0;
			for (size_t i = 0; i < o.rowLength(); i++)
			for (size_t j = 0; j < o.columnLength(); j++)
			storage[cnt++] = o(i, j);
		}
		
		template <class U>
		Matrix(const Matrix<U> &o)
		{
			init(o.rowLength(), o.columnLength());
			size_t cnt = 0;
			for (size_t i = 0; i < o.rowLength(); i++)
			for (size_t j = 0; j < o.columnLength(); j++)
			storage[cnt++] = o(i, j);
		}
		
		Matrix &operator=(const Matrix &o)
		{
			if (this == &o) return *this;
			init(o.rowLength(), o.columnLength());
			size_t cnt = 0;
			for (size_t i = 0; i < o.rowLength(); i++)
			for (size_t j = 0; j < o.columnLength(); j++)
			storage[cnt++] = o(i, j);
			return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
			init(o.rowLength(), o.columnLength());
			size_t cnt = 0;
			for (size_t i = 0; i < o.rowLength(); i++)
			for (size_t j = 0; j < o.columnLength(); j++)
			storage[cnt++] = o(i, j);
			return *this;
		}
		
		Matrix(Matrix &&o) noexcept
		{
			n = std::move(o.n), m = std::move(o.m);
			storage = std::move(o.storage);
			o.storage = nullptr;
		}
		
		Matrix &operator=(Matrix &&o) noexcept
		{
			if (this == &o) return *this;
			if (storage) delete[] storage;
			n = std::move(o.n), m = std::move(o.m);
			storage = std::move(o.storage);
			o.storage = nullptr;
			return *this;
		}
		
		~Matrix() { clear(); }
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
			for (auto iteri = il.begin(); iteri != il.end(); ++iteri)
			if (iteri->size() != il.begin()->size())
			throw std::invalid_argument("invalid initializer_list");
			init(il.size(), il.begin()->size());
			size_t cnt = 0;
			for (auto iteri = il.begin(); iteri != il.end(); ++iteri)
			for (auto iterj = iteri->begin(); iterj != iteri->end(); ++iterj)
			storage[cnt++] = *iterj;
		}
		
		
	private:
		void init(size_t _n, size_t _m, T _init = T())
		{
			n = _n, m = _m;
			size_t tot = n * m;
			if (storage) delete[] storage;
			storage = new T[tot];
			for (size_t i = 0; i < tot; i++) storage[i] = _init;
		}
		
	public:
		size_t rowLength() const { return n; }
		
		size_t columnLength() const { return m; }
		
		void resize(size_t _n, size_t _m, T _init = T())
		{
			if (_n * _m != rowLength() * columnLength())
			{
				T *new_storage = new T[_n * _m];
				size_t tot_new = _n * _m, tot = rowLength() * columnLength();
				for (size_t i = 0; i < tot_new && i < tot; i++)
				new_storage[i] = storage[i];
				for (size_t i = tot; i < tot_new; i++) new_storage[i] = _init;
				delete[] storage;
				storage = new_storage;
			}
			n = _n, m = _m;
		}
		
		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
			resize(sz.first, sz.second, _init);
		}
		
		std::pair<size_t, size_t> size() const
		{
			return std::make_pair(rowLength(), columnLength());
		};
		
		void clear()
		{
			if (storage) delete[] storage;
			n = m = 0;
			storage = nullptr;
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
			if (i >= rowLength() || i < 0 || j >= columnLength() || j < 0)
			throw std::invalid_argument("invalid coordinate");
			return storage[i * m + j];
		}
		
		T &operator()(size_t i, size_t j)
		{
			if (i >= rowLength() || i < 0 || j >= columnLength() || j < 0)
			throw std::invalid_argument("invalid coordinate");
			return storage[i * m + j];
		}
		
		Matrix<T> row(size_t i) const
		{
			if (i >= rowLength() || i < 0)
			throw std::invalid_argument("invalid row");
			
			Matrix<T> result(1, columnLength());
			for (size_t j = 0; j < columnLength(); ++j)
			result(0, j) = storage[i * m + j];
			return result;
		}
		
		Matrix<T> column(size_t i) const
		{
			if (i >= columnLength() || i < 0)
			throw std::invalid_argument("invalid column");
			
			Matrix<T> result(rowLength(), 1);
			for (size_t j = 0; j < rowLength(); ++j)
			result(j, 0) = storage[j * m + i];
			return std::move(result);
		}
		
	public:
		// matlab style
		friend std::ostream &operator<<(std::ostream &out, const Matrix &o)
		{
			int n = o.rowLength(), m = o.columnLength();
			out << '[';
			for (size_t i = 0; i < n - 1; ++i)
			{
				for (size_t j = 0; j < m - 1; ++j) out << o(i, j) << ", ";
				out << o(i, m - 1) << ";\n";
			}
			for (size_t j = 0; j < m - 1; ++j)
			out << o(n - 1, j) << ", ";
			out << o(n - 1, m - 1) << "]";
			return out;
		}
		
	public:
		template <class U>
		bool operator==(const Matrix<U> &o) const
		{
			if (size() != o.size()) return 0;
			for (size_t i = 0; i < rowLength(); i++)
			for (size_t j = 0; j < columnLength(); j++)
			if ((*this)(i, j) != o(i, j)) return 0;
			return 1;
		}
		
		template <class U>
		bool operator!=(const Matrix<U> &o) const
		{
			return !(*this == o);
		}
		
		Matrix operator-() const
		{
			auto result = *this;
			for (size_t i = 0; i < rowLength(); i++)
			for (size_t j = 0; j < columnLength(); j++)
			result(i, j) = -result(i, j);
			return result;
		}
		
		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
			if (size() != o.size())
			throw std::invalid_argument("Sizes of the matrices' don't match");
			for (size_t i = 0; i < rowLength(); i++)
			for (size_t j = 0; j < columnLength(); j++)
			(*this)(i, j) += o(i, j);
			return *this;
		}
		
		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
			*this += (-o);
			return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
			for (size_t i = 0; i < rowLength(); i++)
			for (size_t j = 0; j < columnLength(); j++)
			(*this)(i, j) *= x;
			return *this;
		}
		
		// transpose
		// this function will create a new matrix instead of modifying the original
		Matrix tran() const
		{
			auto result = *this;
			result.resize(columnLength(), rowLength());
			for (size_t i = 0; i < rowLength(); i++)
			for (size_t j = 0; j < columnLength(); j++)
			result(j, i) = (*this)(i, j);
			return result;
		}
		
	private:
		struct IterData
		{
			Matrix *mt{ nullptr };
			std::pair<size_t, size_t> leftUp{ 0, 0 };
			size_t n{ 0 }, m{ 0 }, total{ 0 }, cur{ 0 };
			
			IterData() = default;
			
			IterData(Matrix *_mt, std::pair<size_t, size_t> _lu, size_t _n, size_t _m, size_t _total, size_t _cur) : mt(_mt), leftUp(_lu), n(_n), m(_m), total(_total), cur(_cur) { }
			
			bool operator==(const IterData &o) const
			{
				return mt == o.mt && leftUp == o.leftUp && m == o.m && total == o.total && cur == o.cur;
			}
		};
		
	public: // iterator
		class iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			iterator() = default;
			
			iterator(const iterator &) = default;
			
			iterator &operator=(const iterator &) = default;
			
		private:
			friend class Matrix;
			
			IterData data;
			
			iterator(Matrix *_mt, std::pair<size_t, size_t> _lu, size_t _n, size_t _m, size_t _total, size_t _cur) : data(_mt, _lu, _n, _m, _total, _cur) { }
			
		public:
			difference_type operator-(const iterator &o)
			{
				return data.cur - o.data.cur;
			}
			
			iterator &operator+=(difference_type offset)
			{
				data.cur += offset;
				return *this;
			}
			
			iterator operator+(difference_type offset) const
			{
				auto result = *this;
				return result += offset;
			}
			
			iterator &operator-=(difference_type offset)
			{
				return (*this += (-offset));
			}
			
			iterator operator-(difference_type offset) const
			{
				auto result = *this;
				return (*this -= offset);
			}
			
			iterator &operator++()
			{
				return (*this += 1);
			}
			
			iterator operator++(int)
			{
				auto result = *this;
				++(*this);
				return result;
			}
			
			iterator &operator--()
			{
				return *this -= 1;
			}
			
			iterator operator--(int)
			{
				auto result = *this;
				--(*this);
				return result;
			}
			
			reference operator*() const
			{
				return (*data.mt)(data.leftUp.second + data.cur / data.m, data.leftUp.second + data.cur % data.m);
			}
			
			pointer operator->() const
			{
				return &(operator*());
			}
			
			bool operator==(const iterator &o) const
			{
				return data == o.data;
			}
			
			bool operator!=(const iterator &o) const
			{
				return !(*this == o);
			}
		};
		
		iterator begin()
		{
			return iterator(this, { 0, 0 }, rowLength(), columnLength(), rowLength() * columnLength(), 0);
		}
		
		iterator end()
		{
			return iterator(this, { 0, 0 }, rowLength(), columnLength(), rowLength() * columnLength(), rowLength() * columnLength());
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r)
		{
			auto n = r.first - l.first + 1, m = r.second - l.second + 1;
			return std::make_pair(iterator(this, l, n, m, n * m, 0), iterator(this, l, n, m, n * m, n * m));
		};
		
	public: // policy iterator
		class IterPolicyBase
		{
		protected:
			static void advance(size_t, IterData &);
			
			~IterPolicyBase() = default;
		};
		
		class RowIterator : public IterPolicyBase
		{
		protected:
			static void advance(std::ptrdiff_t n, IterData &data)
			{
				data.cur += n;
			}
		};
		
		class TraceIterator : public IterPolicyBase
		{
		protected:
			static void advance(std::ptrdiff_t n, IterData &data)
			{
				data.cur += n * (1 + data.m);
				if (data.cur == data.total + data.m)
				data.cur = data.total;
			}
		};
		
		class ColumnIterator : public IterPolicyBase
		{
		protected:
			static void advance(std::ptrdiff_t n, IterData &data)
			{
				size_t i = data.cur / data.m;
				size_t j = data.cur % data.m;
				size_t cur = j * data.n + i + n;
				if (cur > data.total)
				data.cur = data.total + 1;
				else if (cur == data.total)
				data.cur = data.total;
				else
				{
					i = cur % data.n;
					j = cur / data.n;
					data.cur = i * data.m + j;
				}
			}
		};
		
		template <class Policy>
		class PolicyIterator : public Policy
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			PolicyIterator() = default;
			
			PolicyIterator(const PolicyIterator &) = default;
			
			PolicyIterator &operator=(const PolicyIterator &) = default;
			
		private:
			friend class Matrix;
			
			IterData data;
			
			PolicyIterator(Matrix *_mt, std::pair<size_t, size_t> _lu, size_t _n, size_t _m, size_t _total, size_t _cur) : data(_mt, _lu, _n, _m, _total, _cur) { }
			
		public:
			difference_type operator-(const PolicyIterator &o)
			{
				return data.cur - o.data.cur;
			}
			
			PolicyIterator &operator+=(difference_type offset)
			{
				Policy::advance(offset, data);
				return *this;
			}
			
			PolicyIterator operator+(difference_type offset) const
			{
				auto result = *this;
				return result += offset;
			}
			
			PolicyIterator &operator-=(difference_type offset)
			{
				return (*this += (-offset));
			}
			
			PolicyIterator operator-(difference_type offset) const
			{
				auto result = *this;
				return (*this -= offset);
			}
			
			PolicyIterator &operator++()
			{
				return (*this += 1);
			}
			
			PolicyIterator operator++(int)
			{
				auto result = *this;
				++(*this);
				return result;
			}
			
			PolicyIterator &operator--()
			{
				return *this -= 1;
			}
			
			PolicyIterator operator--(int)
			{
				auto result = *this;
				--(*this);
				return result;
			}
			
			reference operator*() const
			{
				return (*data.mt)(data.leftUp.second + data.cur / data.m, data.leftUp.second + data.cur % data.m);
			}
			
			pointer operator->() const
			{
				return &(operator*());
			}
			
			bool operator==(const PolicyIterator &o) const
			{
				return data == o.data;
			}
			
			bool operator!=(const PolicyIterator &o) const
			{
				return !(*this == o);
			}
			
			bool operator<(const PolicyIterator &o) const
			{
				return data.cur < o.data.cur;
			}
		};
		
		template <class Policy>
		PolicyIterator<Policy> begin()
		{
			return PolicyIterator<Policy>(this, { 0, 0 }, rowLength(), columnLength(), rowLength() * columnLength(), 0);
		}
		
		template <class Policy>
		PolicyIterator<Policy> end()
		{
			return PolicyIterator<Policy>(this, { 0, 0 }, rowLength(), columnLength(), rowLength() * columnLength(), rowLength() * columnLength());
		}
		
	};
	
}

//
namespace sjtu
{
	
	// scalar multiplication
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)
	{
		Matrix<decltype(T() * U())> result(mat);
		result *= x;
		return result;
	};
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)
	{
		return mat * x;
	};
	
	// multiple
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)
	{
		if (a.columnLength() != b.rowLength())
		throw std::invalid_argument("Matrices are not multipliable");
		
		Matrix<decltype(U() * V())> result(a.rowLength(), b.columnLength());
		for (size_t i = 0; i < a.rowLength(); ++i)
		for (size_t j = 0; j < b.columnLength(); ++j)
		{
			result(i, j) = 0;
			for (size_t k = 0; k < a.columnLength(); ++k)
			result(i, j) += a(i, k) * b(k, j);
		}
		return result;
	};
	
	// add
	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)
	{
		Matrix<decltype(U() * V())> result(a);
		result += b;
		return result;
	};
	
	// sub
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)
	{
		Matrix<decltype(U() * V())> result(a);
		result -= b;
		return result;
	};
	
}

#endif //SJTU_MATRIX_HPP

