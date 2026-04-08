#' Generate block indices for big data grid traversal
#'
#' Partitions a matrix or file into rectangular blocks and returns a data frame
#' of block indices numbered from 1 to K (total number of blocks). Supports
#' exact and flexible partitioning, and multiple traversal orders.
#'
#' @param fnrow Integer. Number of rows in the file or matrix.
#' @param fncol Integer. Number of columns in the file or matrix.
#' @param bnrow Integer. Number of rows per block.
#' @param bncol Integer. Number of columns per block.
#' @param traversal.mode Character. Traversal order for the grid. \code{"rowsnake"}
#'   uses alternating left-right traversal (snake/boustrophedon order);
#'   \code{"lr"} uses strict left-to-right traversal throughout. Defaults to
#'   \code{"rowsnake"}.
#' @param is.flexible Logical or \code{NULL}. If non-\code{NULL}, overrides both
#'   \code{is.flexible.row} and \code{is.flexible.col}. When \code{TRUE}, allows
#'   block dimensions that do not evenly divide the file dimensions. Defaults to
#'   \code{NULL}.
#' @param is.flexible.row Logical. If \code{FALSE}, an error is thrown when
#'   \code{bnrow} does not evenly divide \code{fnrow}. Defaults to \code{FALSE}.
#' @param is.flexible.col Logical. If \code{FALSE}, an error is thrown when
#'   \code{bncol} does not evenly divide \code{fncol}. Defaults to \code{FALSE}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{block_no}{Block index from 1 to K.}
#'     \item{U}{Upper (first) row index of the block.}
#'     \item{D}{Lower (last) row index of the block.}
#'     \item{L}{Left (first) column index of the block.}
#'     \item{R}{Right (last) column index of the block.}
#'   }
#'
#' @examples
#' # Exact partition with snake traversal
#' generate_grid(100, 100, 10, 10)
#'
#' # Flexible partition (block size does not divide file size evenly)
#' generate_grid(105, 103, 10, 10, is.flexible = TRUE)
#'
#' @export
generate_grid <- function(fnrow, fncol, bnrow, bncol, traversal.mode = "rowsnake",
			  is.flexible=NULL, is.flexible.row = FALSE, is.flexible.col = FALSE) {
    if (!is.null(is.flexible)) {
        is.flexible.row <- is.flexible
        is.flexible.col <- is.flexible
    }
    if (!is.flexible.row & (fnrow %% bnrow != 0)) stop(sprintf("Number of rows per block %d must divide number of rows in file %d.", bnrow, fnrow))
    if (!is.flexible.col & (fncol %% bncol != 0)) stop(sprintf("Number of cols per block %d must divide number of cols in file %d.", bncol, fncol))

    if (!(is.flexible.row | is.flexible.col)) return(generate.grid.exact(fnrow, fncol, bnrow, bncol, traversal.mode))

    if (traversal.mode == "rowsnake") return(generate.grid.rowsnake(fnrow, fncol, bnrow, bncol))
    if (traversal.mode == "lr") return(generate.grid.lr(fnrow, fncol, bnrow, bncol))
}

#' Generate an exact block grid analytically
#'
#' Computes block indices when block dimensions divide file dimensions exactly.
#' Uses vectorised operations rather than a loop, making it fast for large grids.
#'
#' @param fnrow Integer. Number of rows in the file or matrix.
#' @param fncol Integer. Number of columns in the file or matrix.
#' @param bnrow Integer. Number of rows per block. Must divide \code{fnrow} exactly.
#' @param bncol Integer. Number of columns per block. Must divide \code{fncol} exactly.
#' @param traversal.mode Character. \code{"rowsnake"} or \code{"lr"}. See
#'   \code{\link{generate_grid}} for details.
#'
#' @return A \code{data.frame} with columns \code{block_no}, \code{U}, \code{D},
#'   \code{L}, \code{R}. See \code{\link{generate_grid}} for column descriptions.
#'
#' @seealso \code{\link{generate_grid}}, \code{\link{generate.grid.rowsnake}},
#'   \code{\link{generate.grid.lr}}
#'
#' @export
generate.grid.exact <- function(fnrow, fncol, bnrow, bncol, traversal.mode) {
    Kr <- fnrow/bnrow
    Kc <- fncol/bncol
    K <- Kr * Kc

    out <- data.frame(block_no = seq(1,K), U = rep(seq(1, fnrow, bnrow), each=Kc), D = rep(seq(bnrow, fnrow, bnrow), each=Kc),
		      L = rep(seq(1, fncol, bncol), Kr), R = rep(seq(bncol, fncol, bncol), Kr))

    if (traversal.mode == "rowsnake") {
        # flip the left-right of every even row (UD), i.e. where block_no is in a certain range.
        out[(out$D %% (2 * bnrow) == 0), "L"] <- rev(rep(seq(1, fncol, bncol), Kr/2))
        out[(out$D %% (2 * bnrow) == 0), "R"] <- rev(rep(seq(bncol, fncol, bncol), Kr/2))
    }

    return(out)
}

#' Generate a flexible block grid with snake traversal
#'
#' Iteratively computes block indices using a snake (boustrophedon) traversal
#' order. Handles cases where block dimensions do not evenly divide file
#' dimensions by allowing the last block in each row or column to be smaller.
#'
#' @param fnrow Integer. Number of rows in the file or matrix.
#' @param fncol Integer. Number of columns in the file or matrix.
#' @param bnrow Integer. Number of rows per block.
#' @param bncol Integer. Number of columns per block.
#'
#' @return A \code{data.frame} with columns \code{block_no}, \code{U}, \code{D},
#'   \code{L}, \code{R}. See \code{\link{generate_grid}} for column descriptions.
#'
#' @seealso \code{\link{generate_grid}}, \code{\link{generate.grid.exact}},
#'   \code{\link{generate.grid.lr}}
#'
#' @export
generate.grid.rowsnake <- function(fnrow, fncol, bnrow, bncol) {
    Kr <- ceiling(fnrow/bnrow)
    Kc <- ceiling(fncol/bncol)

    K <- Kr * Kc

    U <- 1
    D <- U + bnrow - 1
    L <- 1
    R <- L + bncol - 1

    out <- data.frame(block_no = c(1), U = c(U), D = c(D), L = c(L), R = c(R))
    # the direction to traverse
    #TODO: toggle in an option in the future on where to start.
    traversal <- "right"
    traversal.last <- "right"

    #TODO: There should be a way to analytically calculate each partition of the grid without
    #		regenerating the dataframe each time or throwing this in a for loop.
    # 	Check by running testthat.
    for (k in 2:K) {
        if (traversal == "right") {
            L <- min(L + bncol, fncol)
            R <- min(R + bncol, fncol)
            # cannot go further
            if (R == fncol) {
	        traversal <- "down"
                traversal.last <- "right"
	    }
        } else if (traversal == "down") {
            # update columns
            U <- min(U + bnrow, fnrow)
            D <- min(D + bnrow, fnrow)
            if (traversal.last == "right") traversal <- "left"
	    if (traversal.last == "left")  traversal <- "right"

        } else if (traversal == "left") {
            R <- ceiling((R - bncol)/bncol) * bncol
            L <- R - bncol + 1

            # cannot go further
            if (L == 1) {
	        traversal <- "down"
                traversal.last <- "left"
	    }
        }
        out <- rbind(out, data.frame(block_no = c(k), U = c(U), D = c(D), L = c(L), R = c(R)))
    }
    return(out)
}

#' Generate a flexible block grid with left-to-right traversal
#'
#' Iteratively computes block indices using a strict left-to-right traversal
#' order across all rows. Handles cases where block dimensions do not evenly
#' divide file dimensions by allowing the last block in each row or column to
#' be smaller.
#'
#' @param fnrow Integer. Number of rows in the file or matrix.
#' @param fncol Integer. Number of columns in the file or matrix.
#' @param bnrow Integer. Number of rows per block.
#' @param bncol Integer. Number of columns per block.
#'
#' @return A \code{data.frame} with columns \code{block_no}, \code{U}, \code{D},
#'   \code{L}, \code{R}. See \code{\link{generate_grid}} for column descriptions.
#'
#' @seealso \code{\link{generate_grid}}, \code{\link{generate.grid.exact}},
#'   \code{\link{generate.grid.rowsnake}}
#'
#' @export
generate.grid.lr <- function(fnrow, fncol, bnrow, bncol) {
    Kr <- ceiling(fnrow/bnrow)
    Kc <- ceiling(fncol/bncol)

    K <- Kr * Kc

    U <- 1
    D <- U + bnrow - 1
    L <- 1
    R <- L + bncol - 1

    out <- data.frame(block_no = c(1), U = c(U), D = c(D), L = c(L), R = c(R))
    for (k in 2:K) {
	if (R == fncol) {
            U <- min(U + bnrow, fnrow)
            D <- min(D + bnrow, fnrow)
            L <- 1
	    R <- L + bncol - 1
	} else {
            L <- min(L + bncol, fncol)
            R <- min(R + bncol, fncol)
	}
        out <- rbind(out, data.frame(block_no = c(k), U = c(U), D = c(D), L = c(L), R = c(R)))
    }
    return(out)
}
