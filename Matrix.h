#pragma once
#include  <vector>


class Matrix final {
private:
	std::vector<std::vector<double>> matrix;
	size_t cols{};
	size_t rows{};

	typedef unsigned long size_t;
	bool multithreading = false;
	friend class CalculationManager;
	size_t size() const { return cols * rows; };
	const size_t maxRowsSum{ 300 };
	const size_t maxRowsMult{ 200 };
	const size_t maxRowsDet{ 2 };

public:
	Matrix() = default;
	Matrix(const Matrix&) = default;
	Matrix(Matrix&&) = default;
	Matrix(size_t rank);
	Matrix(size_t iRows, size_t iCols);
	Matrix(size_t iRows, size_t iCols, double value);

	Matrix& operator=(const Matrix& iMat) = default;
	Matrix& operator=(Matrix&& iMat) = default;
	~Matrix() = default;

	void setMultithreading(bool status = true);
	void fillWith(double value);
	size_t getRows() const { return rows; }
	Matrix substractWith(const Matrix& iAnother) const;
	Matrix fastSubstractWith(const Matrix& iAnother, size_t numThreads) const;
	Matrix sumWith(const Matrix& iAnother) const;
	Matrix fastSumWith(const Matrix& iAnother, size_t numThreads) const;
	Matrix multiplyWith(const Matrix& another) const;
	Matrix fastMultiplyWith(const Matrix& iAnother, size_t numThreads) const;

	double det() const;
	double simpleDet(const Matrix& iMat) const;
	double fastDet(const Matrix& iMat, size_t numThreads) const;
	Matrix minor(const Matrix& iMat, size_t rowIndex, size_t colIndex) const;

	static Matrix createDiagonal(size_t rank, double value= 1.0);

	friend Matrix operator+(const Matrix& iFirst, const Matrix& iSecond);
	friend Matrix operator-(const Matrix& iFirst, const Matrix& iSecond);
	friend Matrix operator*(const Matrix& iFirst, const Matrix& iSecond);

	bool operator==(const Matrix& another) const
	{
		if (rows != another.rows || cols != another.cols) return false;
		for (size_t i = 0; i < rows; i++)
		{
			if (matrix[i] != another.matrix[i]) return false;
		}
		return true;
	}
};

class CalculationManager
{
public:

	CalculationManager(const Matrix& iM1, const Matrix& iM2, size_t numThreads) :numThreads(numThreads), m1(iM1), m2(iM2) {};
	CalculationManager(const Matrix& iM1, size_t numThreads) :numThreads(numThreads), m1(iM1), m2(iM1) {};
	~CalculationManager() = default;
	Matrix sum();
	Matrix multiply();
	Matrix substract();
	double det();
private:
	size_t numThreads;
	const Matrix& m1;
	const Matrix& m2;
	std::vector<std::pair<size_t, size_t>> makeIntervalsForDet() const;
	std::vector<std::pair<size_t, size_t>> makeIntervals() const;

	Matrix calculate(void (CalculationManager::* f)(std::pair<size_t, size_t>&, Matrix*));
	void subSum(std::pair<size_t, size_t>& interval, Matrix* ioResult);
	void subMulti(std::pair<size_t, size_t>& interval, Matrix* ioResult);
	void subSubstr(std::pair<size_t, size_t>& interval, Matrix* ioResult);
	double subDet(std::pair<size_t, size_t>& interval);
	double calculateDet();
};