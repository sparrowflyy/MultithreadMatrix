#include "pch.h"
#include "Matrix.h"
#include <chrono>

static const unsigned long big_size_for_mult = 2300;
static const unsigned long big_size_for_sum = 7000;

static Matrix big_mult_matrix(big_size_for_mult, big_size_for_mult, 666);
static Matrix big_diagonal_matrix = Matrix::createDiagonal(big_size_for_mult, 1);
static Matrix big_filled_matrix(big_size_for_sum, big_size_for_sum, 666);
static Matrix big_empty_matrix(big_size_for_sum, big_size_for_sum);

void timeTestSum(Matrix& A, Matrix& B, size_t thread_count)
{
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	Matrix C = A.fastSumWith(B, thread_count);
	std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Parallel Execution of sum " << A.getRows() << "x" << A.getRows() << " with " << thread_count << " threads : " << dur.count() << "s" << std::endl;
}

void timeTestMult(Matrix& A, Matrix& B, size_t thread_count)
{
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	Matrix C = A.fastMultiplyWith(B, thread_count);
	std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Parallel Execution of mult " << A.getRows() << "x" << A.getRows() << " with " << thread_count << " threads : " << dur.count() <<"s"<< std::endl;
}
void timeTestDet(Matrix& A, size_t thread_count)
{
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	double det = A.fastDet(A, thread_count);
	std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Parallel Execution of det " << A.getRows() << "x" << A.getRows() << " with " << thread_count << " threads : " << dur.count() << "s" << std::endl;
}

TEST(Matrix_determinant, simple) {
	Matrix identity = Matrix::createDiagonal(11, 2);											   
	double det = identity.det();																   
	EXPECT_TRUE(det == 2048);
}
TEST(Matrix_determinant, multithreading) {
	Matrix identity = Matrix::createDiagonal(11, 2);
	identity.setMultithreading();
	double det = identity.det();
	EXPECT_TRUE(det == 2048);
}

TEST(Matrix_multiplication, simple) {

	Matrix multiplied = big_diagonal_matrix * big_mult_matrix;
	EXPECT_TRUE(multiplied == big_mult_matrix);
}
TEST(Matrix_multiplication, multithreading) {
	big_diagonal_matrix.setMultithreading();
	Matrix multiplied = big_diagonal_matrix * big_mult_matrix;
	big_diagonal_matrix.setMultithreading(false);
	EXPECT_TRUE(multiplied == big_mult_matrix);
}
TEST(Matrix_sum, simple) {

	Matrix sum = big_empty_matrix + big_filled_matrix;
	EXPECT_TRUE(sum == big_filled_matrix);
}
TEST(Matrix_sum, multithreading) {
	big_empty_matrix.setMultithreading();
	Matrix sum = big_empty_matrix + big_filled_matrix;
	big_empty_matrix.setMultithreading(false);
	EXPECT_TRUE(sum == big_filled_matrix);
}
TEST(Matrix_substract, simple) {
	Matrix substract = big_filled_matrix - big_empty_matrix;
	EXPECT_TRUE(substract == big_filled_matrix);
}
TEST(Matrix_substract, multithreading) {
	big_filled_matrix.setMultithreading();
	Matrix substract = big_filled_matrix - big_empty_matrix;
	EXPECT_TRUE(substract == big_filled_matrix);
}
TEST(time_research, sum5000)
{
	Matrix A(5000, 5000, 1);
	Matrix B(5000, 5000, 1);
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	Matrix C = A.sumWith(B);
	std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
	std::cout << "Execution of sum " << A.getRows() << "x" << A.getRows() << " " << dur.count() << std::endl;

	for (size_t i = 2; i <= 10; i++)
	{
		timeTestSum(A, B, i);
	}
}

TEST(time_research, mult2000)
{
	Matrix A(2000, 2000, 1);
	Matrix B(2000, 2000, 1);
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	Matrix C = A.multiplyWith(B);
	std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
	std::cout << "Execution of sum " << A.getRows() << "x" << A.getRows() << ": " << dur.count() <<"s"<<std::endl;

	for (size_t i = 2; i <= 10; i++)
	{
		timeTestMult(A, B, i);
	}
}
TEST(time_research, det10x10)
{
	Matrix A(10, 10, 1);

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	double det = A.det();
	std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
	std::cout << "Execution of det " << A.getRows() << "x" << A.getRows() << " " << dur.count() << "s" << std::endl;

	for (size_t i = 2; i <= 8; i+=2)
	{
		timeTestDet(A, i);
	}
}