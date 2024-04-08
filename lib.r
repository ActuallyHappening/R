pow = function(x, n) Reduce(`%*%`, replicate(n, x, simplify = FALSE))
`%^%` = pow

# prints all its arguments
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
	if (all(A %*% t(A) == diag(n))) {
		return(TRUE)
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
	if (all(A %*% A == A)) {
		return(TRUE)
	}
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
				print_all("Matrix converges to a matrix with rows 1=2")
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
	} else {
		print("Conclusion: Matrix is not Markov (doesn't sum to 1)")
	}
}