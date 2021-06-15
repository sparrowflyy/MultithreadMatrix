#include  "Matrix.h"
#include <thread>
#include <future>
#include <iostream>
#include <deque>


Matrix::Matrix(size_t iRows, size_t iCols)
{
	rows = iRows;
	cols = iCols;
	matrix.resize(rows);
	for (size_t i = 0; i < rows; i++)
	{
		matrix[i].resize(cols);
		std::fill(matrix[i].begin(), matrix[i].end(), 0.0);
	}
}

Matrix::Matrix(size_t iRows, size_t iCols, double value)
{
	rows = iRows;
	cols = iCols;
	matrix.resize(rows);
	for (size_t i = 0; i < rows; i++)
	{
		matrix[i].resize(cols);
		std::fill(matrix[i].begin(), matrix[i].end(), 0.0);
	}
}

void Matrix::setMultithreading(bool status)
{
	multithreading = status;
}

Matrix operator+(const Matrix& iFirst, const Matrix& iSecond)
{
	if (!iFirst.multithreading  && !iSecond.multithreading ) {
		return iFirst.sumWith(iSecond);
	}
	return iFirst.fastSumWith(iSecond, std::thread::hardware_concurrency());
}

Matrix operator*(const Matrix& iFirst, const Matrix& iSecond)
{
	if (!iFirst.multithreading && !iSecond.multithreading) {
		return iFirst.multiplyWith(iSecond);
	}
	return iFirst.fastMultiplyWith(iSecond,std::thread::hardware_concurrency());
}

Matrix operator-(const Matrix& iFirst, const Matrix& iSecond)
{
	if (!iFirst.multithreading && !iSecond.multithreading) {
		return iFirst.substractWith(iSecond);
	}
	return iFirst.fastSubstractWith(iSecond, std::thread::hardware_concurrency());
}

double Matrix::det() const
{
	if (multithreading) 
		return fastDet(*this, std::thread::hardware_concurrency());
	return simpleDet(*this);
}


Matrix Matrix::createDiagonal(size_t rank, double value)
{
	Matrix diagonal(rank, rank);
	for (size_t i = 0; i < rank; i++) {
		diagonal.matrix[i][i] = value;
	}
	return diagonal;
}

void Matrix::fillWith(double value)
{
	for (size_t i = 0; i < rows; i++) {
		std::fill(matrix[i].begin(), matrix[i].end(), value);
	}
}
Matrix Matrix::sumWith(const Matrix& iAnother) const
{
	if (rows != iAnother.rows || cols != iAnother.cols)
		return Matrix();
	Matrix sum(rows, cols);
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			sum.matrix[i][j] = matrix[i][j] + iAnother.matrix[i][j];
		}
	}
	return sum;
}

Matrix Matrix::substractWith(const Matrix& iAnother) const
{
	if (rows != iAnother.rows || cols != iAnother.cols) return Matrix();
	Matrix substruct(rows, cols);
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			substruct.matrix[i][j] = matrix[i][j] - iAnother.matrix[i][j];
		}
	}
	return substruct;
}

Matrix Matrix::multiplyWith(const Matrix& iAnother) const
{
	Matrix result(rows, iAnother.cols);
	for (size_t i = 0; i < result.rows; i++)
	{
		for (size_t j = 0; j < result.cols; j++)
		{
			for (size_t k = 0; k < result.cols; k++)
			{
				result.matrix[i][j] += matrix[i][k] * iAnother.matrix[k][j];
			}
		}
	}
	return result;
}
Matrix::Matrix(size_t rank)
{
	rows = rank;
	cols = rank;
	matrix.resize(rows);
	for (size_t i = 0; i < rows; i++)
	{
		matrix[i].resize(cols);
	}
}


Matrix Matrix::minor(const Matrix& iMat, size_t rowIndex, size_t colIndex) const
{
	Matrix sub_mat(iMat.rows - 1);
	for (size_t i = 1; i < iMat.rows; i++)
	{
		const std::vector<double>& temp_row = iMat.matrix[i];
		std::copy(temp_row.begin(), temp_row.begin() + colIndex, sub_mat.matrix[i - 1].begin());
		std::copy(temp_row.begin() + colIndex + 1, temp_row.end(), sub_mat.matrix[i - 1].begin() + colIndex);

	}
	return sub_mat;
}

double Matrix::simpleDet(const Matrix& mat) const
{
	double det = 0;
	int sign = 1;
	if (mat.cols == 1) return mat.matrix[0][0];
	if (mat.cols == 2) return mat.matrix[0][0] * mat.matrix[1][1] - mat.matrix[0][1] * mat.matrix[1][0];
	for (size_t i = 0; i < mat.cols; i++) {
		det = det + (sign * mat.matrix[0][i] * simpleDet(minor(mat, 0, i)));
		sign = -sign;
	}
	return det;
}


double CalculationManager::subDet(std::pair<size_t, size_t>& iInterval)
{
	double det = 0;
	int sign;
	if (iInterval.first % 2 == 0) sign = 1;
	else sign = -1;
	for (size_t i = iInterval.first; i < iInterval.second; i++) {
		det = det + (sign * m1.matrix[0][i] * m1.simpleDet(m1.minor(m1, 0, i)));
		sign = -sign;
	}
	return det;
}

void CalculationManager::subSubstr(std::pair<size_t, size_t>& interval, Matrix* ioResult)
{
#ifdef DEBUG
	std::cout << "Thread (" << interval.first << ";" << interval.second << ") starts work\n";
#endif // DEBUG

	for (size_t i = interval.first; i < interval.second; i++)
	{
		for (size_t j = 0; j < ioResult->cols; j++)
		{

			ioResult->matrix[i][j] = m1.matrix[i][j] - m2.matrix[i][j];
		}
	}

#ifdef DEBUG
	std::cout << "Thread (" << interval.first << ";" << interval.second << ") end work\n";
#endif // DEBUG
	
}

void CalculationManager::subSum(std::pair<size_t, size_t>& interval, class Matrix* ioResult)
{
	for (size_t i = interval.first; i < interval.second; i++)
	{
		for (size_t j = 0; j < ioResult->cols; j++)
		{

			ioResult->matrix[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
		}
	}
}
void CalculationManager::subMulti(std::pair<size_t, size_t>& iInterval, Matrix* ioResult)
{
	for (size_t i = iInterval.first; i < iInterval.second; i++)
	{
		for (size_t j = 0; j < ioResult->cols; j++)
		{
			for (size_t k = 0; k < ioResult->cols; k++)
			{
				ioResult->matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
			}
		}
	}
}

std::vector<std::pair<size_t, size_t>> CalculationManager::makeIntervals() const
{
	std::vector<std::pair<size_t, size_t>> intervals;
	size_t distance = m1.rows / numThreads;
	size_t current_row = 0;
	for (size_t i = 0; i < numThreads; i++)
	{
		intervals.push_back({ current_row,current_row + distance });
		current_row += distance;

	}
	intervals[numThreads - 1].second = m1.rows;
	return intervals;
}
std::vector<std::pair<size_t, size_t>> CalculationManager::makeIntervalsForDet() const
{
	std::vector<std::pair<size_t, size_t>> intervals;
	size_t distance = m1.rows / numThreads;
	if (distance == 0 || m1.rows % numThreads != 0) distance++;
	size_t current_row = 0;
	for (; current_row + distance < m1.rows; current_row += distance) {
		intervals.push_back({ current_row,current_row + distance });
	}
	intervals.push_back({ current_row,m1.rows });
	return intervals;
}


double CalculationManager::calculateDet()
{
	auto intervals = makeIntervalsForDet();
	std::vector<std::future<double>> all_futures;
	for (size_t i = 0; i < intervals.size(); i++) {
		all_futures.push_back(std::async(std::launch::async, &CalculationManager::subDet, this, std::ref(intervals[i])));
	}
	double det = 0;
	for (size_t i = 0; i < intervals.size(); i++) {
		det += all_futures[i].get();
	}
	return det;
}

Matrix CalculationManager::calculate(void (CalculationManager::* f)(std::pair<size_t, size_t>&, Matrix*))
{
	Matrix result(m1.rows, m1.cols);
	auto intervals = makeIntervals();
	std::vector<std::future<void>> all_futures(intervals.size() - 1);
	for (size_t i = 0; i < all_futures.size(); i++) {
		std::packaged_task<void(std::pair<size_t, size_t>&, Matrix*)> temp_task([&](std::pair<size_t, size_t>& a, Matrix* b) { return (this->*f)(a, b); });
		all_futures[i] = std::async(std::launch::async, std::move(temp_task), std::ref(intervals[i]), &result);
	}
	(this->*f)(intervals[intervals.size() - 1], &result);
	for (size_t i = 0; i < intervals.size() - 1; i++) {
		all_futures[i].get();
	}
	return result;

}
double CalculationManager::det()
{
	return calculateDet();
}
Matrix CalculationManager::multiply()
{
	return calculate(&CalculationManager::subMulti);
}
Matrix CalculationManager::sum()
{
	return calculate(&CalculationManager::subSum);
}
Matrix CalculationManager::substract()
{
	return calculate(&CalculationManager::subSubstr);
}

Matrix Matrix::fastMultiplyWith(const Matrix& another, size_t num_of_threads) const
{
	if (cols != another.rows) return Matrix();
	size_t threads_count = rows / maxRowsMult + 1;
	if (threads_count == 1) return multiplyWith(another);
	if (threads_count > num_of_threads) threads_count = num_of_threads;
	CalculationManager multiplier(*this, another, threads_count);
	return multiplier.multiply();
}
double Matrix::fastDet(const Matrix& mat, size_t num_of_threads) const
{
	size_t threads_count = rows / maxRowsDet;
	if (threads_count == 1) return simpleDet(mat);
	if (threads_count > num_of_threads) threads_count = num_of_threads;
	CalculationManager detCalculator(mat, threads_count);
	return detCalculator.det();
}
Matrix Matrix::fastSubstractWith(const Matrix& another, size_t num_of_threads) const
{
	if (rows != another.rows || cols != another.cols) return Matrix();
	size_t threads_count = rows / maxRowsSum + 1;
	if (threads_count == 1) return substractWith(another);
	if (threads_count > num_of_threads) threads_count = num_of_threads;
	CalculationManager Substructor(*this, another, threads_count);
	return Substructor.substract();
}

Matrix Matrix::fastSumWith(const Matrix& another, size_t num_of_threads) const
{
	if (rows != another.rows || cols != another.cols) return Matrix();
	size_t threads_count = rows / maxRowsSum + 1;
	if (threads_count == 1) return sumWith(another);
	if (threads_count > num_of_threads) threads_count = num_of_threads;
	CalculationManager adder(*this, another, threads_count);
	return adder.sum();
}