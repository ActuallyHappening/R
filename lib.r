pow = function(x, n) Reduce(`%*%`, replicate(n, x, simplify = FALSE))
`%^%` = pow

#' prints all its arguments
print_all = function(...) {
	for (i in 1:length(list(...))) {
		print(list(...)[[i]])
	}
}

## Is orthonormal rotation matrix

is_orthonormal <- function(A) {
	n <- nrow(A)
	if (n != ncol(A)) {
		# must be square matrix
		return(FALSE)
	}
	## A %*% t(A) == I
	I = diag(n)
	LHS = round(A %*% t(A), 4)
	if (all(LHS == I)) {
		print_all("Satisfies first orthonormal condition", "~LHS", LHS, "I", I)
		if (abs(det(A)) == 1) {
			print_all("Satisfies second orthonormal condition", "det(A)", det(A))
			return(TRUE)
		} else {
			print_all("Conclusion: Not an orthonormal rotation matrix: det(A) != 1", "det(A)", det(A))
		}
	} else {
		print_all("Conclusion: Not an orthonormal projection matrix: A %*% t(A) != I", "~LHS", LHS, "I", I)
	}
	return(FALSE)
}

## is perpendicular projection matrix

is_perpendicular_projection <- function(A) {
	n <- nrow(A)
	if (n != ncol(A)) {
		# must be square matrix
		return(FALSE)
	}
	A2 = A %*% A
	if (all(A2 == A)) {
		print_all("Conclusion: Matrix is perpendicular projection")
		return(TRUE)
	}
	print_all("Conclusion: Matrix is not a perpendicular projection", "A2", A2, "A", A)
	return(FALSE)
}

## is markov transition matrix, showing that all elements are non-negative and row sum is 1
## and prints out convergence (or lack of convergence) for powers 128, 256, and 512, 513, and 514

is_markov = function(A) {
	n = nrow(A)
	if (n != ncol(A)) {
		# must be square matrix
		return(FALSE)
	}

	if (all(A >= 0) & all(A <= 1) & all(rowSums(A) == 1)) {
		print("Matrix is Markov")

		A512 = A %^% 512
		A513 = A512 %*% A
		A514 = A513 %*% A
		A515 = A514 %*% A
		print_all("A512", A512, "A513", A513, "A514", A514, "A515", A515)

		if (all(A512 == A513) & all(A513 == A514)) {
			print("Matrix converges to a specific matrix")
			# check if rows have different values
			if (all(A512[1,] == A512[2,])) {
				print_all("Matrix converges to a matrix with rows 1=2, hence doesn't depend on starting state")
			} else {
				print_all("Matrix converges, but it depends what your starting state was")
			}
		} else {
			print_all("Conclusion: Matrix does not converge to a single value")
			if (all(A512 == A514)) {
				print_all("Conclusion: Matrix oscillates with a frequency of 2")
			} else if (all(A512 == A515)) {
				print_all("Conclusion: Matrix oscillates with a frequency of 3")
			} else {
				print_all("Conclusion: Matrix does not converge with a nice oscillation")
			}
		}

		return(TRUE)
	} else {
		print("Conclusion: Matrix is not Markov (doesn't sum to 1)")
	}
	return(FALSE)
}

analyze_matrix = function(M) {
	print_all("Analyzing matrix", M)

	markov = is_markov(M)
	orthonormal = is_orthonormal(M)
	projection = is_perpendicular_projection(M)

	if (markov) {
		print("Summary: Matrix is Markov")
	}
	if (orthonormal) {
		print("Summary: Matrix is orthonormal")
	}
	if (projection) {
		print("Summary: Matrix is a perpendicular projection")
	}

	if (!markov & !orthonormal & !projection) {
		print("Summary: Matrix is not Markov, orthonormal, or a perpendicular projection")
	}
}