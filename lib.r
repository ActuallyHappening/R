# install.packages("diagram")
# mirror 5

# all ( c(TRUE, TRUE) ) == TRUE
# all ( c(TRUE, FALSE) ) == FALSE
# all ( matrix(c(TRUE, TRUE, TRUE, TRUE), nrow=2, byrow=TRUE) ) == TRUE
# nrow(M) == 69, ncol(M) = 420


# rownames(M) = c("A", "B")
# colnames(M) = rownames(M)
# plotmat(M, pos = nrow(M))

# π(0)P5 =(1,0,0,0,0)P5

# rep(1,5) == c(1,1,1,1,1)
# replicate(5, M, simplify = FALSE) is kind of list(M, M, M, M, M)

pow = function(x, n) Reduce(`%*%`, replicate(n, x, simplify = FALSE))
`%^%` = pow

recursive_pow = function(A, n) {
	if (n == 0) {
		return(diag(nrow(A)))
	}
	if (n == 1) {
		return(A)
	}
	return(A %*% recursive_pow(A, n-1))
}

loop_pow = function(A, n) {
	result = diag(nrow(A))
	for (i in 1:n) {
		result = result %*% A
	}
	return(result)
}

#' prints all its arguments
print_all = function(...) {
	for (i in 1:length(list(...))) {
		print(list(...)[[i]])
	}
}

## Is orthonormal rotation matrix

#region is_orthonormal
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
#endregion

## is perpendicular projection matrix

#region is_perpendicular_projection
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
#endregion

## is markov transition matrix, showing that all elements are non-negative and row sum is 1
## and prints out convergence (or lack of convergence) for powers 128, 256, and 512, 513, and 514

#region is_markov
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
#endregion

#region analyze_matrix
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
#endregion

extrapolate_markov = function (M, i) {
	print_all("Analyzing markov matrix", M, "with initial starting vector / state", i)

	for (n in 2:10) {
		exp = 2 ** n
		Mn = M %^% exp

		v = i %*% Mn
		print_all("Extrapolating to power", exp, "gives Mn", Mn, "with final value v =", v)
	}
}

analyze_markov = function(M) {
	print("Is markov?")
	is_markov(M)
	print("Assuming is markov.")

	n = nrow(M)
	for (start_i in 1:n) {
		# create starting vector with all 0s except for 1 at start_i
		i = rep(0, n)
		i[start_i] = 1

		extrapolate_markov(M, i)
	}
}

# initial<-matrix(c(1,0,0,0,0),ncol=5) #start in Bronze P<-matrix(c(0.6, 0.25, 0.1, 0.04, 0.01,
# P<-matrix(c(0.6, 0.25, 0.1, 0.04, 0.01,
#   0.4, 0.25, 0.3, 0.04, 0.01,
#   0, 0.4, 0.4, 0.18, 0.02,
#   0, 0, 0.4, 0.56, 0.04,
#   0, 0, 0, 0.8, 0.2),
# ncol=5,byrow=TRUE) #transition matrix